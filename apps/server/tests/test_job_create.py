"""Jobs from FASTQ files in the data folder."""

import time
from pathlib import Path
from typing import Any

import pytest
from fastapi.testclient import TestClient
from scalehd.simulate_run import RunSample, write_run
from scalehd_server.app import create_app
from scalehd_server.config import ServerSettings

MODEL = {"method": "model"}


def _job(*samples: dict[str, Any], name: str = "run 01") -> dict[str, Any]:
    return {"name": name, "samples": list(samples), "settings": MODEL}


def test_a_job_from_the_data_folder_is_queued(
    queued: TestClient, data_root: Path, pools: Any
) -> None:
    body = _job(
        {"name": "a", "r1": "run-01/a_R1.fastq.gz", "r2": "run-01/a_R2.fastq.gz"},
        {"name": "b", "r1": "run-01/b_R1.fastq.gz"},
    )
    response = queued.post("/api/jobs", json=body)
    assert response.status_code == 201
    job = response.json()
    assert job["name"] == "run 01"
    assert job["demo"] is False
    assert job["method"] == "model"
    run = data_root.resolve() / "run-01"
    assert [(s["name"], s["r1"], s["r2"]) for s in job["samples"]] == [
        ("a", str(run / "a_R1.fastq.gz"), str(run / "a_R2.fastq.gz")),
        ("b", str(run / "b_R1.fastq.gz"), None),
    ]
    assert (Path(job["output_dir"]) / "job.json").is_file()
    assert len(pools.submitted) == 1  # one worker for test
    assert queued.get("/api/jobs").json()[0]["id"] == job["id"]


@pytest.mark.parametrize(
    ("sample", "problem"),
    [
        ({"name": "a", "r1": "run-01/missing_R1.fastq.gz"}, "no file"),
        ({"name": "a", "r1": "run-01/a_R1.fastq.gz", "r2": "run-01/missing.fastq.gz"}, "no file"),
        ({"name": "a", "r1": "../outside.fastq.gz"}, "no file"),
        ({"name": "a", "r1": "run-01"}, "no file"),  # a folder, not a file
    ],
)
def test_files_must_be_in_the_data_folder(
    queued: TestClient, sample: dict[str, Any], problem: str
) -> None:
    response = queued.post("/api/jobs", json=_job(sample))
    assert response.status_code == 400
    assert response.json()["detail"].startswith("sample a: ")
    assert problem in response.json()["detail"]
    assert queued.get("/api/jobs").json() == []


def test_a_server_without_a_data_folder_refuses_jobs(logged_in: TestClient) -> None:
    response = logged_in.post("/api/jobs", json=_job({"name": "a", "r1": "a_R1.fastq.gz"}))
    assert response.status_code == 409
    assert "data folder" in response.json()["detail"]


def test_a_simulated_samples_truth_file_is_used(queued: TestClient, data_root: Path) -> None:
    write_run(
        data_root / "sim", [RunSample("s", ("17_1_1_7_2", "43_1_1_7_2"))], pairs=10, workers=1
    )
    body = _job({"name": "s", "r1": "sim/s_S1_L001_R1_001.fastq.gz"})
    job = queued.post("/api/jobs", json=body).json()
    assert job["samples"][0]["truth"] == "17_1_1_7_2/43_1_1_7_2"


def test_a_real_job_runs_end_to_end(tmp_path: Path) -> None:
    root = tmp_path / "data"
    samples = [
        RunSample("paired-17-43", ("17_1_1_7_2", "43_1_1_7_2")),
        RunSample("r1-only-20-21", ("20_1_1_7_2", "21_1_1_7_2"), single_end=True),
    ]
    write_run(root / "run", samples, pairs=1500, seed=4, workers=1)
    settings = ServerSettings(
        database_dir=tmp_path / "db", workspace=tmp_path / "ws", data_root=root, workers=1
    )
    with TestClient(create_app(settings)) as client:
        body = {"username": "autotest-user", "password": "correct horse"}
        assert client.post("/api/auth/register", json=body).status_code == 201
        listing = client.get("/api/inputs", params={"folder": "run"}).json()
        chosen = [
            {"name": s["name"], "r1": s["r1"], "r2": s["r2"]}
            for s in listing["samples"]
            if not s["undetermined"] and not s["skipped"]
        ]
        job = client.post("/api/jobs", json=_job(*chosen)).json()
        deadline = time.monotonic() + 180
        while job["status"] != "finished":
            assert time.monotonic() < deadline, f"job still {job['status']}: {job['samples']}"
            time.sleep(0.25)
            job = client.get(f"/api/jobs/{job['id']}").json()

    assert [s["name"] for s in job["samples"]] == ["paired-17-43", "r1-only-20-21"]
    for sample in job["samples"]:
        assert sample["status"] == "finished", sample["error"]
        assert sample["genotype"] == sample["truth"]
        assert sample["matches_truth"] is True
        assert (Path(job["output_dir"]) / sample["name"] / "call.json").is_file()
