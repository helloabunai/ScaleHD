"""Job and sample rows as API bodies."""

from pathlib import Path

from scalehd_server.models import Job, Sample, SampleStatus, User
from scalehd_server.progress import ProgressReporter
from scalehd_server.schemas import GenotypeMethod, JobSettings, job_out, job_summary
from sqlalchemy.orm import Session, sessionmaker


def test_job_out_counts_progress_and_carries_the_truth(sessions: sessionmaker[Session]) -> None:
    with sessions() as session:
        job = Job(
            owner=User(username="autotest-user", password_hash="x"),
            name="Demo: 9 simulated samples",
            demo=True,
            settings=JobSettings(method=GenotypeMethod.MODEL).model_dump(mode="json"),
            output_dir="/ws/autotest-user/1-demo",
        )
        job.samples = [
            Sample(
                name="normal-heterozygote",
                status=SampleStatus.FINISHED,
                genotype="17_1_1_7_2/21_1_1_7_2",
                confidence=40.0,
                flags=["low_depth"],
                simulation={"alleles": ["17_1_1_7_2", "21_1_1_7_2"], "pairs": 5000, "seed": 1},
                truth="17_1_1_7_2/21_1_1_7_2",
                matches_truth=True,
            ),
            Sample(name="expanded", status=SampleStatus.FAILED, error="ValueError: x"),
            Sample(name="homozygous"),
        ]
        session.add(job)
        session.commit()
        out = job_out(job)
        summary = job_summary(job)

    assert (out.sample_count, out.samples_done) == (3, 2)
    assert out.demo
    assert out.method is GenotypeMethod.MODEL
    assert out.output_dir == summary.output_dir == "/ws/autotest-user/1-demo"
    first = out.samples[0]
    assert first.r1 is None
    assert (first.truth, first.matches_truth) == ("17_1_1_7_2/21_1_1_7_2", True)
    assert out.samples[2].status is SampleStatus.QUEUED
    assert out.model_dump(mode="json")["samples"][1]["error"] == "ValueError: x"


def test_job_out_says_how_far_running_samples_have_got(
    sessions: sessionmaker[Session], tmp_path: Path
) -> None:
    folder = tmp_path / "1-run"
    with sessions() as session:
        job = Job(
            owner=User(username="autotest-user", password_hash="x"),
            name="run",
            settings=JobSettings().model_dump(mode="json"),
            output_dir=str(folder),
        )
        job.samples = [
            Sample(name="counting", status=SampleStatus.RUNNING),
            Sample(name="starting", status=SampleStatus.RUNNING),
            Sample(name="done", status=SampleStatus.FINISHED),
        ]
        session.add(job)
        session.commit()
        for name in ("counting", "done"):
            ProgressReporter(folder / name).counted(0.25)
        out = job_out(job)

    samples = out.model_dump(mode="json")["samples"]
    assert samples[0]["progress"] == {"stage": "counting", "share": 0.25}
    assert samples[1]["progress"] is None
    assert samples[2]["progress"] is None
