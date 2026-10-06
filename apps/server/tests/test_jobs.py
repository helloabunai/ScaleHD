"""The demo job and the jobs API, plus one real run end to end."""

import json
import os
import shutil
import time
from pathlib import Path
from typing import Any

import pytest
from fastapi.testclient import TestClient
from scalehd_server import demo
from scalehd_server.app import create_app
from scalehd_server.config import ServerSettings
from scalehd_server.demo import DEMO_SAMPLES, DemoSample
from scalehd_server.models import Job, JobStatus
from sqlalchemy import update


def test_demo_job_queues_thirteen_simulated_samples(quiet_client: TestClient, pools: Any) -> None:
    response = quiet_client.post("/api/jobs/demo")
    assert response.status_code == 201
    job = response.json()
    assert job["name"] == "Demo: 13 simulated samples"
    assert job["demo"] is True
    assert job["method"] == "model"
    assert job["sample_count"] == 13
    assert [s["name"] for s in job["samples"]] == [s.name for s in DEMO_SAMPLES]
    truths = {s["name"]: s["truth"] for s in job["samples"]}
    assert truths["loss-of-interruption"] == "19_1_1_7_2/42_0_1_7_2"
    assert truths["homozygous"] == "21_1_1_7_2/21_1_1_7_2"
    assert truths["ccg-7-and-10"] == "17_1_1_7_2/17_1_1_10_2"
    assert truths["caacag-duplication"] == "19_2_1_10_2/40_1_1_7_2"
    assert truths["ccgcca-deletion"] == "19_1_0_7_2/40_1_1_7_2"
    assert truths["ccgcca-duplication"] == "17_1_1_7_2/42_1_2_7_2"
    assert truths["ccgcca-deletion-and-insertion"] == "19_1_0_7_2/42_1_2_7_2"
    assert len(pools.submitted) == 1  # one worker: one sample handed out


def test_demo_job_writes_its_folder(quiet_client: TestClient, tmp_path: Path) -> None:
    job = quiet_client.post("/api/jobs/demo").json()
    folder = (tmp_path / "workspace" / "autotest-user" / f"{job['id']}-demo").resolve()
    assert job["output_dir"] == str(folder)
    record = json.loads((folder / "job.json").read_text())
    assert record["owner"] == "autotest-user"
    assert record["settings"]["method"] == "model"
    assert len(record["samples"]) == 13


def test_demo_uses_the_model_method_even_when_the_default_is_legacy(
    quiet_client: TestClient,
) -> None:
    quiet_client.put("/api/settings", json={"method": "legacy"})
    assert quiet_client.get("/api/settings").json()["method"] == "legacy"
    assert quiet_client.post("/api/jobs/demo").json()["method"] == "model"


def test_running_the_demo_twice_makes_two_jobs(quiet_client: TestClient) -> None:
    first = quiet_client.post("/api/jobs/demo").json()
    second = quiet_client.post("/api/jobs/demo").json()
    assert first["id"] != second["id"]
    assert first["output_dir"] != second["output_dir"]
    listed = quiet_client.get("/api/jobs").json()
    assert [job["id"] for job in listed] == [second["id"], first["id"]]


def test_jobs_are_private(quiet_client: TestClient) -> None:
    job = quiet_client.post("/api/jobs/demo").json()
    bob = TestClient(quiet_client.app)
    bob.post("/api/auth/register", json={"username": "bob", "password": "correct horse"})
    assert bob.get("/api/jobs").json() == []
    someone_elses = bob.get(f"/api/jobs/{job['id']}")
    missing = bob.get("/api/jobs/999999")
    assert someone_elses.status_code == missing.status_code == 404
    assert someone_elses.json() == missing.json() == {"detail": "no such job"}


def test_jobs_need_a_login(client: TestClient) -> None:
    assert client.post("/api/jobs/demo").status_code == 401
    assert client.get("/api/jobs/1").status_code == 401
    assert client.post("/api/jobs/1/cancel").status_code == 401
    assert client.delete("/api/jobs/1").status_code == 401


@pytest.mark.skipif(os.geteuid() == 0, reason="root can write anywhere")
def test_an_unwritable_workspace_fails_with_the_folder_named(
    quiet_client: TestClient, tmp_path: Path
) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir(exist_ok=True)
    workspace.chmod(0o500)
    try:
        response = quiet_client.post("/api/jobs/demo")
    finally:
        workspace.chmod(0o700)
    assert response.status_code == 500
    assert response.json()["detail"].startswith(
        f"can't write to the workspace at {workspace.resolve()}/autotest-user/"
    )
    assert quiet_client.get("/api/jobs").json() == []


def test_the_demo_runs_end_to_end(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        demo,
        "DEMO_SAMPLES",
        (
            DemoSample("normal-heterozygote", ("17_1_1_7_2", "21_1_1_7_2"), seed=1, pairs=1000),
            DemoSample("expanded", ("17_1_1_7_2", "43_1_1_7_2"), seed=2, pairs=1000),
        ),
    )
    # One worker for two samples: the second is handed out when the first finishes.
    settings = ServerSettings(
        database_dir=tmp_path / "data", workspace=tmp_path / "workspace", workers=1
    )
    with TestClient(create_app(settings)) as client:
        client.post(
            "/api/auth/register", json={"username": "autotest-user", "password": "correct horse"}
        )
        job = client.post("/api/jobs/demo").json()
        assert job["name"] == "Demo: 2 simulated samples"
        deadline = time.monotonic() + 120
        while job["status"] != "finished":
            assert time.monotonic() < deadline, f"demo still {job['status']}: {job['samples']}"
            time.sleep(0.25)
            job = client.get(f"/api/jobs/{job['id']}").json()

    assert job["samples_done"] == 2
    for sample in job["samples"]:
        assert sample["status"] == "finished", sample["error"]
        assert sample["genotype"] == sample["truth"]
        assert sample["matches_truth"] is True
        assert sample["confidence"] > 0
        folder = Path(job["output_dir"]) / sample["name"]
        assert (folder / "counts.json").exists()
        assert (folder / "call.json").exists()
        assert (folder / "input" / f"{sample['name']}_R1.fastq.gz").exists()


def _set_job(client: TestClient, job_id: int, **values: Any) -> None:
    """Change a stored job directly, e.g. mark it finished without running it."""
    with client.app.state.sessions() as session:
        session.execute(update(Job).where(Job.id == job_id).values(**values))
        session.commit()


def test_deleting_a_finished_job_removes_its_folder_and_records(
    quiet_client: TestClient,
) -> None:
    job = quiet_client.post("/api/jobs/demo").json()
    _set_job(quiet_client, job["id"], status=JobStatus.FINISHED)
    folder = Path(job["output_dir"])
    assert folder.exists()
    assert quiet_client.delete(f"/api/jobs/{job['id']}").status_code == 204
    assert not folder.exists()
    assert quiet_client.get(f"/api/jobs/{job['id']}").status_code == 404
    assert quiet_client.get("/api/jobs").json() == []


def test_a_running_job_cannot_be_deleted(quiet_client: TestClient) -> None:
    job = quiet_client.post("/api/jobs/demo").json()
    response = quiet_client.delete(f"/api/jobs/{job['id']}")
    assert response.status_code == 409
    assert response.json() == {"detail": "the job is still running"}
    assert Path(job["output_dir"]).exists()


def test_someone_elses_job_cannot_be_deleted(quiet_client: TestClient) -> None:
    job = quiet_client.post("/api/jobs/demo").json()
    _set_job(quiet_client, job["id"], status=JobStatus.FINISHED)
    bob = TestClient(quiet_client.app)
    bob.post("/api/auth/register", json={"username": "bob", "password": "correct horse"})
    assert bob.delete(f"/api/jobs/{job['id']}").status_code == 404
    assert Path(job["output_dir"]).exists()


def test_a_job_whose_folder_is_already_gone_can_be_deleted(quiet_client: TestClient) -> None:
    job = quiet_client.post("/api/jobs/demo").json()
    _set_job(quiet_client, job["id"], status=JobStatus.FINISHED)
    shutil.rmtree(job["output_dir"])
    assert quiet_client.delete(f"/api/jobs/{job['id']}").status_code == 204


def test_a_folder_outside_your_workspace_is_never_deleted(
    quiet_client: TestClient, tmp_path: Path
) -> None:
    precious = tmp_path / "precious"
    precious.mkdir()
    job = quiet_client.post("/api/jobs/demo").json()
    _set_job(quiet_client, job["id"], status=JobStatus.FINISHED, output_dir=str(precious))
    response = quiet_client.delete(f"/api/jobs/{job['id']}")
    assert response.status_code == 500
    assert response.json()["detail"].startswith(f"won't delete {precious}")
    assert precious.exists()
    assert [j["id"] for j in quiet_client.get("/api/jobs").json()] == [job["id"]]


def _wait_for(client: TestClient, job_id: int, status: str) -> dict[str, Any]:
    """helper for job status changes"""
    deadline = time.monotonic() + 5
    while (job := client.get(f"/api/jobs/{job_id}").json())["status"] != status:
        assert time.monotonic() < deadline, f"job still {job['status']}"
        time.sleep(0.01)
    return job


def test_a_cancelled_job_can_be_deleted_once_its_running_samples_finish(
    quiet_client: TestClient, pools: Any
) -> None:
    job = quiet_client.post("/api/jobs/demo").json()
    response = quiet_client.post(f"/api/jobs/{job['id']}/cancel")
    assert response.status_code == 200
    cancelling = response.json()
    # One worker/one sample was running, wait for it to finish. The other 12 never start.
    assert cancelling["status"] == "cancelling"
    assert [s["status"] for s in cancelling["samples"]] == ["running"] + ["cancelled"] * 12
    assert quiet_client.post(f"/api/jobs/{job['id']}/cancel").json()["status"] == "cancelling"
    # The running sample still writes into the job folder even under cancel context
    refused = quiet_client.delete(f"/api/jobs/{job['id']}")
    assert refused.status_code == 409
    assert refused.json() == {
        "detail": "the job is waiting for already-processing samples to finish"
    }

    _, future = pools.submitted[0]
    future.set_exception(ValueError("no molecules"))
    cancelled = _wait_for(quiet_client, job["id"], "cancelled")
    assert cancelled["finished_at"] is not None
    assert len(pools.submitted) == 1
    assert quiet_client.delete(f"/api/jobs/{job['id']}").status_code == 204
    assert not Path(job["output_dir"]).exists()


def test_a_job_that_has_stopped_cannot_be_cancelled(quiet_client: TestClient) -> None:
    job = quiet_client.post("/api/jobs/demo").json()
    _set_job(quiet_client, job["id"], status=JobStatus.FINISHED)
    response = quiet_client.post(f"/api/jobs/{job['id']}/cancel")
    assert response.status_code == 409
    assert response.json() == {"detail": "the job has already stopped"}


def test_someone_elses_job_cannot_be_cancelled(quiet_client: TestClient) -> None:
    job = quiet_client.post("/api/jobs/demo").json()
    bob = TestClient(quiet_client.app)
    bob.post("/api/auth/register", json={"username": "bob", "password": "correct horse"})
    assert bob.post(f"/api/jobs/{job['id']}/cancel").status_code == 404
    assert quiet_client.get(f"/api/jobs/{job['id']}").json()["status"] == "running"
