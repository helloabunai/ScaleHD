"""Legacy (ScaleHD 1.x) jobs, their settings, running a sample through the runner, sharing
cores, and reading the results.

The runner here is a stand-in. The real one only runs in the Docker image.
"""

import json
import shutil
import sys
import time
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from fastapi.testclient import TestClient
from pydantic import ValidationError
from scalehd_server import legacy, worker
from scalehd_server.app import create_app
from scalehd_server.config import ServerSettings
from scalehd_server.models import Job, Sample, User
from scalehd_server.progress import SampleStage, read_progress
from scalehd_server.runner import JobRunner
from scalehd_server.schemas import (
    REFUSED_ERROR_TOLERANCES,
    GenotypeMethod,
    JobSettings,
    LegacyDemultiplexFlags,
    LegacySettings,
    LegacyTrimFlags,
)
from scalehd_server.worker import SampleTask
from scalehd_server.workspace import sample_folder_names
from sqlalchemy.orm import Session, sessionmaker

DATA = Path(__file__).parent / "data"
LEGACY_JOB = JobSettings(method=GenotypeMethod.LEGACY).for_job()


def result(name: str = "atypical") -> dict[str, Any]:
    """What the runner saved for a sample in tools/legacy-parity.sh -- "atypical" (a typical
    and a realigned atypical allele) or "few-reads" (1.x failed in atypical scanning)."""
    found: dict[str, Any] = json.loads((DATA / f"legacy-result-{name}.json").read_text())
    return found


# Settings


def test_settings_are_written_as_1x_reads_them() -> None:
    config = LegacySettings().config()

    assert config["instance_flags"]["demultiplex"] == "False"
    assert config["instance_flags"]["quality_control"] == "True"
    assert config["alignment_flags"]["indel_penalty"] == "6,6"
    assert config["alignment_flags"]["chain_drop"] == "0.5"
    assert config["trim_flags"]["error_tolerance"] == "0.39"
    assert config["demultiplex_flags"]["min_length"] == ""
    assert config["demultiplex_flags"]["max_length"] == ""


def test_the_error_tolerances_refused_are_those_1x_refuses() -> None:
    # 1.x's check (its __backend.py), on every hundredth from 0 to 1.
    allowed = np.arange(0, 1.1, 0.01)
    refused = {i / 100 for i in range(101) if float(f"{i / 100:.2f}") not in allowed}

    assert refused == REFUSED_ERROR_TOLERANCES


@pytest.mark.parametrize("tolerance", [*sorted(REFUSED_ERROR_TOLERANCES), 0.355, 1.2])
def test_an_error_tolerance_1x_would_refuse_is_refused(tolerance: float) -> None:
    with pytest.raises(ValidationError):
        LegacyTrimFlags(error_tolerance=tolerance)


def test_other_error_tolerances_are_kept() -> None:
    assert LegacyTrimFlags(error_tolerance=0.36).error_tolerance == 0.36


def test_adapters_are_tidied_and_checked() -> None:
    assert LegacyDemultiplexFlags(forward_adapter=" gcgaccctgg ").forward_adapter == "GCGACCCTGG"
    with pytest.raises(ValidationError):
        LegacyDemultiplexFlags(forward_adapter="GCGNCC")


def test_only_the_adapter_positions_1x_handles_are_offered() -> None:
    # AP is in 1.x's check, but batchadapt then gets no flag and writes nothing.
    with pytest.raises(ValidationError):
        LegacyDemultiplexFlags(forward_position="AP")  # type: ignore[arg-type]


@pytest.mark.parametrize("stage", ["sequence_alignment", "genotype_prediction"])
def test_stages_1x_cant_do_without_stay_on(stage: str) -> None:
    with pytest.raises(ValidationError):
        LegacySettings.model_validate({"instance_flags": {stage: False}})


def test_a_job_keeps_1xs_settings_only_when_it_runs_1x() -> None:
    with_both = JobSettings(legacy=LegacySettings(cores_per_sample=4))

    assert with_both.for_job().legacy is None
    legacy_job = with_both.model_copy(update={"method": GenotypeMethod.LEGACY}).for_job()
    assert legacy_job.legacy is not None
    assert legacy_job.legacy.cores_per_sample == 4
    assert JobSettings(method=GenotypeMethod.LEGACY).for_job().legacy == LegacySettings()


# The call


def test_the_call_is_1xs_report_row() -> None:
    found = result()
    row = found["instance_report"].splitlines()[1].split(",")

    call = legacy.call(found)

    assert call["schema"] == "scalehd.legacy-call/1"
    assert call["genotype"] == f"{row[2]}/{row[17]}" == "17_1_1_7_2/22_2_1_7_2"
    for allele, start in zip(call["alleles"], (2, 17), strict=True):
        columns = row[start : start + 14]
        assert allele["structure"] == columns[0]
        assert allele["status"] == columns[1]
        assert [str(allele[key]) for key in ("map_fw", "map_percent_fw", "purged_fw")] == (
            columns[2:5]
        )
        assert [str(allele[key]) for key in ("map_rv", "map_percent_rv", "purged_rv")] == (
            columns[5:8]
        )
        assert str(allele["backward_slippage"]) == columns[8]
        assert str(allele["somatic_mosaicism"]) == columns[9]
        assert allele["variant_call"] == columns[10]
        assert allele["intervening"] == columns[12]
        assert str(allele["confidence"]) == columns[13]
    assert [str(value) for value in call["flag_values"].values()] == row[33:50]


def test_an_allele_has_its_structure_counts() -> None:
    atypical = legacy.call(result())["alleles"][1]

    assert atypical["status"] == "Atypical"
    assert (atypical["cag"], atypical["caacag"], atypical["ccgcca"]) == (22, 2, 1)
    assert (atypical["ccg"], atypical["cct"]) == (7, 2)


def test_raised_flags_leave_out_the_one_that_says_filtering_passed() -> None:
    call = legacy.call(result())

    assert call["flag_values"][legacy.PASS_FLAG] is True
    assert call["flags"] == [
        column
        for column, raised in call["flag_values"].items()
        if raised and column != legacy.PASS_FLAG
    ]


def test_the_samples_confidence_is_its_least_confident_alleles() -> None:
    call = legacy.call(result())
    call["alleles"][0]["confidence"] = 72

    assert legacy.confidence(call) == 72.0


def test_matching_the_truth_ignores_allele_order() -> None:
    call = legacy.call(result())

    assert legacy.matches(call, "22_2_1_7_2/17_1_1_7_2") is True
    assert legacy.matches(call, "17_1_1_7_2/22_1_1_7_2") is False


def test_cag_reads_come_from_the_alleles_own_distribution() -> None:
    sample = result()["sample"]
    typical, atypical = sample["primaryallele"], sample["secondaryallele"]

    for allele, cag in ((typical, 17), (atypical, 22)):
        reads = legacy.cag_reads(allele)
        assert len(reads) == 200
        assert reads.index(max(reads)) + 1 == cag
    # Realigned: 1.x's own reference first, then the typical one.
    assert len(atypical["fwarray"]) == 8000


def test_a_failure_in_1xs_words() -> None:
    assert legacy.failure(result("few-reads")) == (
        "Legacy ScaleHD failed in atypical scanning: "
        "<200 aligned reads in Allele #1. Data un-usable."
    )


# Running a sample


def fake_runner(tmp_path: Path, saved: str | None = "atypical", code: int = 0) -> Path:
    """Stands in for the legacy runner: reports a few stages, then saves a result."""
    runner = tmp_path / "fake-runner"
    source = str(DATA / f"legacy-result-{saved}.json") if saved else None
    runner.write_text(
        f"""#!{sys.executable}
import json, shutil, sys
from pathlib import Path

out = Path(sys.argv[2])
if out.exists() and any(out.iterdir()):
    sys.exit(2)
for stage in ("preparing", "alignment", "genotyping"):
    print(json.dumps({{"stage": stage}}), flush=True)
print("1.x was here", file=sys.stderr)
if {source!r}:
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy({source!r}, out / "result.json")
sys.exit({code})
"""
    )
    runner.chmod(0o755)
    return runner


def _task(tmp_path: Path, runner: Path, **change: Any) -> SampleTask:
    reads = tmp_path / "reads"
    reads.mkdir(exist_ok=True)
    for name in ("x_R1.fastq.gz", "x_R2.fastq.gz"):
        (reads / name).write_bytes(b"")
    task = {
        "sample_id": 1,
        "name": "sample 1",
        "folder": tmp_path / "job" / "sample-1",
        "settings": LEGACY_JOB,
        "r1": reads / "x_R1.fastq.gz",
        "r2": reads / "x_R2.fastq.gz",
        "threads": 3,
        "legacy_runner": runner,
        "truth": "17_1_1_7_2/22_2_1_7_2",
    }
    return SampleTask(**(task | change))


def test_a_legacy_sample_is_run_through_the_runner(tmp_path: Path) -> None:
    task = _task(tmp_path, fake_runner(tmp_path))

    found = worker.run_task(task)

    assert found.genotype == "17_1_1_7_2/22_2_1_7_2"
    assert found.confidence == 100.0
    assert found.matches_truth is True
    assert found.call is not None
    assert found.call["schema"] == "scalehd.legacy-call/1"
    assert found.flags == found.call["flags"]
    request = json.loads((task.folder / legacy.REQUEST).read_text())
    # "sample 1" has a space 1.x can't take, so it's run under its folder's name.
    assert request["label"] == "sample-1"
    assert request["threads"] == 3
    assert request["config"] == LegacySettings().config()
    assert "1.x was here" in (task.folder / legacy.LOG).read_text()
    progress = read_progress(task.folder)
    assert progress is not None
    assert progress.stage is SampleStage.GENOTYPING
    assert progress.share == legacy.stage_progress(LegacySettings(), 3)[SampleStage.GENOTYPING]


def test_an_earlier_attempts_output_is_cleared_first(tmp_path: Path) -> None:
    task = _task(tmp_path, fake_runner(tmp_path))
    (task.folder / legacy.FOLDER).mkdir(parents=True)
    (task.folder / legacy.FOLDER / "old").write_text("")

    assert worker.run_task(task).genotype == "17_1_1_7_2/22_2_1_7_2"


def test_a_sample_1x_failed_says_why(tmp_path: Path) -> None:
    task = _task(tmp_path, fake_runner(tmp_path, "few-reads", code=1))

    with pytest.raises(legacy.LegacyError, match="failed in atypical scanning: <200 aligned"):
        worker.run_task(task)


def test_a_runner_that_saved_nothing(tmp_path: Path) -> None:
    task = _task(tmp_path, fake_runner(tmp_path, None, code=2))

    with pytest.raises(legacy.LegacyError, match=r"stopped \(exit 2\): 1.x was here"):
        worker.run_task(task)


def test_1x_needs_both_reads(tmp_path: Path) -> None:
    task = _task(tmp_path, fake_runner(tmp_path), r2=None)

    with pytest.raises(legacy.LegacyError, match="both read files"):
        worker.run_task(task)


def test_progress_fills_by_a_fixed_part_per_stage() -> None:
    done = legacy.stage_progress(LegacySettings())

    order = [
        SampleStage.PREPARING,
        SampleStage.INDEXING,
        SampleStage.QUALITY_CONTROL,
        SampleStage.ALIGNMENT,
        SampleStage.ATYPICAL_SCANNING,
        SampleStage.GENOTYPING,
        SampleStage.SNP_CALLING,
    ]
    shares = [done[stage] for stage in order]
    assert shares[0] == 0
    assert shares == sorted(shares)
    assert shares[-1] < 1
    # Realignment isn't counted ahead: progress holds while it runs.
    assert done[SampleStage.REALIGNMENT] == done[SampleStage.GENOTYPING]
    assert SampleStage.DEMULTIPLEXING not in done


def test_progress_counts_only_the_stages_the_job_runs() -> None:
    settings = LegacySettings.model_validate(
        {"instance_flags": {"demultiplex": True, "quality_control": False, "snp_calling": False}}
    )

    done = legacy.stage_progress(settings)

    assert SampleStage.DEMULTIPLEXING in done
    assert SampleStage.QUALITY_CONTROL not in done
    assert SampleStage.SNP_CALLING not in done
    # Genotyping is then the last stage, so has the rest to itself.
    seconds = {
        stage: legacy.stage_seconds(stage, 1)
        for stage in done
        if stage is not SampleStage.REALIGNMENT
    }
    genotyping = seconds[SampleStage.GENOTYPING]
    assert done[SampleStage.GENOTYPING] == pytest.approx(1 - genotyping / sum(seconds.values()))


def test_progress_follows_the_cores_a_sample_gets() -> None:
    alignment, scanning = SampleStage.ALIGNMENT, SampleStage.ATYPICAL_SCANNING

    assert legacy.stage_seconds(alignment, 4) < legacy.stage_seconds(alignment, 1) / 3
    assert legacy.stage_seconds(scanning, 2) == legacy.stage_seconds(scanning, 1)
    assert legacy.stage_seconds(scanning, 8) == legacy.stage_seconds(scanning, 3)
    assert legacy.stage_seconds(scanning, 3) < legacy.stage_seconds(scanning, 1)
    one, many = (
        legacy.stage_progress(LegacySettings(), 1),
        legacy.stage_progress(LegacySettings(), 8),
    )
    assert one[SampleStage.GENOTYPING] > many[SampleStage.GENOTYPING]


def test_names_1x_can_use_are_kept(tmp_path: Path) -> None:
    assert legacy.label("HD_123-a", tmp_path / "hd-123-a") == "HD_123-a"
    assert legacy.label("HD.123", tmp_path / "hd-123") == "hd-123"


# Sharing cores


def _jobs(
    sessions: sessionmaker[Session], tmp_path: Path, *settings: JobSettings, samples: int = 2
) -> list[int]:
    with sessions() as session:
        owner = User(username="autotest-user", password_hash="x")
        jobs = []
        for n, job_settings in enumerate(settings):
            job = Job(
                owner=owner,
                name=f"job {n}",
                settings=job_settings.model_dump(mode="json"),
                output_dir=str(tmp_path / "ws" / f"{n}-job"),
            )
            job.samples = [
                Sample(name=f"s{i}", r1=f"/r/s{i}_R1.fq", r2=f"/r/s{i}_R2.fq")
                for i in range(samples)
            ]
            jobs.append(job)
        session.add_all(jobs)
        session.commit()
        return [job.id for job in jobs]


def _legacy(cores: int) -> JobSettings:
    return JobSettings(method=GenotypeMethod.LEGACY, legacy=LegacySettings(cores_per_sample=cores))


def test_legacy_samples_take_their_cores_and_the_queue_waits_its_turn(
    sessions: sessionmaker[Session], tmp_path: Path, pools: Any
) -> None:
    legacy_job, model_job = _jobs(sessions, tmp_path, _legacy(3), JobSettings())
    runner = JobRunner(sessions, workers=4, executor_factory=pools, legacy_runner=Path("/r"))
    runner.submit(legacy_job)
    runner.submit(model_job)

    # 3 of 4 cores taken: the next legacy sample needs 3, and the model ones wait behind it.
    assert [(task.name, task.threads) for task, _ in pools.submitted] == [("s0", 3)]
    first, future = pools.submitted[0]
    assert first.legacy_runner == Path("/r")
    future.set_exception(legacy.LegacyError("no"))

    def started() -> list[int]:
        return [task.threads for task, _ in pools.submitted]

    deadline = time.monotonic() + 5
    while len(started()) < 3:
        assert time.monotonic() < deadline, "timed out waiting"
        time.sleep(0.01)
    assert started() == [3, 3, 1]


def test_a_sample_never_waits_for_more_cores_than_there_are(
    sessions: sessionmaker[Session], tmp_path: Path, pools: Any
) -> None:
    (job,) = _jobs(sessions, tmp_path, _legacy(16))
    JobRunner(sessions, workers=4, executor_factory=pools).submit(job)

    assert [task.threads for task, _ in pools.submitted] == [4]


def test_model_samples_take_one_core_each(
    sessions: sessionmaker[Session], tmp_path: Path, pools: Any
) -> None:
    (job,) = _jobs(sessions, tmp_path, JobSettings(), samples=3)
    JobRunner(sessions, workers=2, executor_factory=pools).submit(job)

    assert [task.threads for task, _ in pools.submitted] == [1, 1]


# Through the API


@pytest.fixture
def legacy_server(tmp_path: Path, pools: Any) -> Iterator[TestClient]:
    """Logged in, with a data folder of two paired samples and one R1 only, and a
    (stand-in) legacy runner."""
    data = tmp_path / "data"
    (data / "run").mkdir(parents=True)
    for name in (
        "a_R1.fastq.gz",
        "a_R2.fastq.gz",
        "b_R1.fastq.gz",
        "b_R2.fastq.gz",
        "c_R1.fastq.gz",
    ):
        (data / "run" / name).write_bytes(b"")
    settings = ServerSettings(
        database_dir=tmp_path / "db",
        workspace=tmp_path / "workspace",
        data_root=data,
        workers=4,
        legacy_runner=fake_runner(tmp_path),
    )
    with TestClient(create_app(settings, executor_factory=pools)) as client:
        body = {"username": "autotest-user", "password": "correct horse"}
        assert client.post("/api/auth/register", json=body).status_code == 201
        yield client


def _pair(name: str) -> dict[str, str]:
    return {"name": name, "r1": f"run/{name}_R1.fastq.gz", "r2": f"run/{name}_R2.fastq.gz"}


def _legacy_job(client: TestClient, **settings: Any) -> dict[str, Any]:
    body = {
        "name": "1.x run",
        "samples": [_pair("a"), _pair("b")],
        "settings": {"method": "legacy", **settings},
    }
    response = client.post("/api/jobs", json=body)
    assert response.status_code == 201, response.text
    job: dict[str, Any] = response.json()
    return job


def test_the_server_offers_legacy_when_it_has_the_runner(legacy_server: TestClient) -> None:
    options = legacy_server.get("/api/settings/options").json()

    assert (options["methods"], options["cores"]) == (["legacy", "model"], 4)


def test_a_legacy_job_keeps_1xs_settings(legacy_server: TestClient, pools: Any) -> None:
    job = _legacy_job(legacy_server, legacy={"cores_per_sample": 2})

    assert job["method"] == "legacy"
    assert job["settings"]["legacy"]["cores_per_sample"] == 2
    assert job["settings"]["legacy"]["trim_flags"]["error_tolerance"] == 0.39
    assert [task.threads for task, _ in pools.submitted] == [2, 2]


def test_a_model_job_doesnt_keep_1xs_settings(legacy_server: TestClient) -> None:
    body = {"name": "new", "samples": [_pair("a")], "settings": {"method": "model"}}

    assert legacy_server.post("/api/jobs", json=body).json()["settings"]["legacy"] is None


def test_a_legacy_job_needs_both_reads_of_each_sample(legacy_server: TestClient) -> None:
    body = {
        "name": "1.x run",
        "samples": [_pair("a"), {"name": "c", "r1": "run/c_R1.fastq.gz"}],
        "settings": {"method": "legacy"},
    }
    response = legacy_server.post("/api/jobs", json=body)

    assert response.status_code == 400
    assert response.json()["detail"].endswith("these have one: c")


def test_bad_1x_settings_are_refused(legacy_server: TestClient) -> None:
    body = {
        "name": "1.x run",
        "samples": [_pair("a")],
        "settings": {"method": "legacy", "legacy": {"trim_flags": {"error_tolerance": 0.35}}},
    }

    assert legacy_server.post("/api/jobs", json=body).status_code == 422


def _ran(job: dict[str, Any], *saved: str) -> list[Path]:
    """Put a saved result in each sample's folder, as if the runner had run them."""
    names = sample_folder_names([sample["name"] for sample in job["samples"]])
    folders = [Path(job["output_dir"]) / name for name in names]
    for folder, name in zip(folders, saved, strict=True):
        (folder / legacy.FOLDER).mkdir(parents=True)
        shutil.copy(DATA / f"legacy-result-{name}.json", folder / legacy.FOLDER / legacy.RESULT)
    return folders


def test_1xs_report_for_the_whole_job(legacy_server: TestClient) -> None:
    job = _legacy_job(legacy_server)
    _ran(job, "atypical", "few-reads")

    response = legacy_server.get(f"/api/jobs/{job['id']}/legacy/InstanceReport.csv")

    assert response.status_code == 200
    assert response.headers["content-disposition"].endswith('-InstanceReport.csv"')
    header, *rows = response.text.splitlines()
    assert header.startswith("SampleName,,Primary GTYPE,Status")
    assert [row.split(",")[0] for row in rows] == ["atypical", "few-reads"]
    distributions = legacy_server.get(f"/api/jobs/{job['id']}/legacy/AlignedDistributions.csv")
    # The failed sample has no distributions.
    assert {row.split(",")[0] for row in distributions.text.splitlines()[1:]} == {"atypical"}


def test_no_1x_report_before_anything_has_run(legacy_server: TestClient) -> None:
    job = _legacy_job(legacy_server)

    assert legacy_server.get(f"/api/jobs/{job['id']}/legacy/InstanceReport.csv").status_code == 404


def test_a_legacy_samples_page(legacy_server: TestClient) -> None:
    job = _legacy_job(legacy_server)
    _ran(job, "atypical", "few-reads")
    first, failed = (sample["id"] for sample in job["samples"])

    detail = legacy_server.get(f"/api/jobs/{job['id']}/samples/{first}").json()["legacy"]
    assert detail["status"] == "complete"
    assert detail["call"]["genotype"] == "17_1_1_7_2/22_2_1_7_2"
    # Both alleles at CCG 7, but the atypical one was realigned: a chart each.
    charts = detail["cag_charts"]
    assert [(chart["alleles"], chart["called"], chart["ccg"]) for chart in charts] == [
        (["17_1_1_7_2"], [17], 7),
        (["22_2_1_7_2"], [22], 7),
    ]
    assert all(len(chart["reads"]) == 200 for chart in charts)
    assert all(len(chart["reads"]) == 20 for chart in detail["ccg_charts"])

    detail = legacy_server.get(f"/api/jobs/{job['id']}/samples/{failed}").json()["legacy"]
    assert detail["status"] == "failed"
    assert detail["failed_stage"] == "DSP"
    assert detail["call"] is None


def test_a_legacy_samples_files(legacy_server: TestClient) -> None:
    job = _legacy_job(legacy_server)
    (folder, _) = _ran(job, "atypical", "few-reads")
    sample = job["samples"][0]["id"]
    detail = legacy_server.get(f"/api/jobs/{job['id']}/samples/{sample}").json()

    assert detail["files"] == ["r1", "r2", "legacy-result"]
    response = legacy_server.get(f"/api/jobs/{job['id']}/samples/{sample}/files/legacy-result")
    assert response.status_code == 200
    assert response.headers["content-disposition"].endswith(f'{folder.name}-legacy-result.json"')
    # 1.x's graphs are outside this sample's folder in the saved result, so not offered.
    url = f"/api/jobs/{job['id']}/samples/{sample}/files/legacy-graphs"
    assert legacy_server.get(url).status_code == 404


def test_a_legacy_jobs_settings_in_the_validation_export(legacy_server: TestClient) -> None:
    job = _legacy_job(legacy_server)
    text = legacy_server.get(f"/api/jobs/{job['id']}/validation.csv").text

    assert "legacy.trim_flags.error_tolerance=0.39" in text
    assert "legacy.alignment_flags.indel_penalty=[6, 6]" in text
    assert "min_molecules" not in text
