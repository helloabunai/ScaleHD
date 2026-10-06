"""Job scheduling, recording, failures and restarts, with stand-in pools."""

import json
import subprocess
import sys
import textwrap
import time
from collections.abc import Callable
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path
from typing import Any

from scalehd_server.models import Job, JobStatus, Sample, SampleStatus, User
from scalehd_server.runner import JobRunner
from scalehd_server.schemas import GenotypeMethod, JobSettings
from scalehd_server.worker import SampleResult
from sqlalchemy.orm import Session, sessionmaker

Q, R, F, X = SampleStatus.QUEUED, SampleStatus.RUNNING, SampleStatus.FINISHED, SampleStatus.FAILED
C = SampleStatus.CANCELLED


def make_job(
    sessions: sessionmaker[Session], tmp_path: Path, n: int, status: SampleStatus = Q
) -> int:
    with sessions() as session:
        job = Job(
            owner=User(username="autotest-user", password_hash="x"),
            name="Demo",
            demo=True,
            status=JobStatus.RUNNING if status == R else JobStatus.QUEUED,
            settings=JobSettings(method=GenotypeMethod.MODEL).model_dump(mode="json"),
            output_dir=str(tmp_path / "ws" / "autotest-user" / "1-demo"),
        )
        job.samples = [
            Sample(
                name=f"s{i}",
                status=status,
                simulation={"alleles": ["21_1_1_7_2"], "pairs": 100, "seed": i},
                truth="21_1_1_7_2/21_1_1_7_2",
            )
            for i in range(n)
        ]
        session.add(job)
        session.commit()
        return job.id


def statuses(sessions: sessionmaker[Session], job_id: int) -> tuple[JobStatus, list[SampleStatus]]:
    with sessions() as session:
        job = session.get(Job, job_id)
        assert job is not None
        return job.status, [sample.status for sample in job.samples]


def eventually(check: Callable[[], bool], timeout: float = 5) -> None:
    """Results are recorded on their own thread, so wait for the expected state."""
    deadline = time.monotonic() + timeout
    while not check():
        assert time.monotonic() < deadline, "timed out waiting"
        time.sleep(0.01)


def result_for(task: Any) -> SampleResult:
    return SampleResult(
        sample_id=task.sample_id,
        call={"genotype": "21_1_1_7_2/21_1_1_7_2"},
        genotype="21_1_1_7_2/21_1_1_7_2",
        confidence=42.0,
        flags=["homozygous"],
        matches_truth=True,
    )


def test_one_worker_runs_one_sample_at_a_time(
    sessions: sessionmaker[Session], tmp_path: Path, pools: Any
) -> None:
    job_id = make_job(sessions, tmp_path, 3)
    JobRunner(sessions, workers=1, executor_factory=pools).submit(job_id)
    assert len(pools.submitted) == 1
    assert statuses(sessions, job_id) == (JobStatus.RUNNING, [R, Q, Q])

    task, future = pools.submitted[0]
    future.set_result(result_for(task))
    eventually(lambda: len(pools.submitted) == 2)
    assert statuses(sessions, job_id) == (JobStatus.RUNNING, [F, R, Q])

    task, future = pools.submitted[1]
    future.set_result(result_for(task))
    eventually(lambda: len(pools.submitted) == 3)
    task, future = pools.submitted[2]
    future.set_result(result_for(task))
    eventually(lambda: statuses(sessions, job_id)[0] == JobStatus.FINISHED)
    assert statuses(sessions, job_id) == (JobStatus.FINISHED, [F, F, F])
    with sessions() as session:
        job = session.get(Job, job_id)
        assert job is not None
        assert job.started_at is not None
        assert job.finished_at is not None
        first = job.samples[0]
        assert (first.genotype, first.confidence, first.flags) == (
            "21_1_1_7_2/21_1_1_7_2",
            42.0,
            ["homozygous"],
        )
        assert first.matches_truth is True
        assert first.call == {"genotype": "21_1_1_7_2/21_1_1_7_2"}


def test_a_task_describes_its_sample(
    sessions: sessionmaker[Session], tmp_path: Path, pools: Any
) -> None:
    job_id = make_job(sessions, tmp_path, 1)
    JobRunner(sessions, workers=1, executor_factory=pools).submit(job_id)
    task, _ = pools.submitted[0]
    assert task.folder == tmp_path / "ws" / "autotest-user" / "1-demo" / "s0"
    assert task.settings.method is GenotypeMethod.MODEL
    assert task.simulation == {"alleles": ["21_1_1_7_2"], "pairs": 100, "seed": 0}
    assert task.truth == "21_1_1_7_2/21_1_1_7_2"
    assert task.r1 is None


def test_a_failed_sample_is_recorded_and_the_job_still_finishes(
    sessions: sessionmaker[Session], tmp_path: Path, pools: Any
) -> None:
    job_id = make_job(sessions, tmp_path, 2)
    JobRunner(sessions, workers=2, executor_factory=pools).submit(job_id)
    (first, first_future), (second, second_future) = pools.submitted
    first_future.set_exception(ValueError("no molecules"))
    second_future.set_result(result_for(second))
    eventually(lambda: statuses(sessions, job_id)[0] == JobStatus.FINISHED)
    assert statuses(sessions, job_id) == (JobStatus.FINISHED, [X, F])
    with sessions() as session:
        sample = session.get(Sample, first.sample_id)
        assert sample is not None
        assert sample.error == "ValueError: no molecules"


def test_a_crashed_worker_fails_its_sample_and_a_fresh_pool_is_made(
    sessions: sessionmaker[Session], tmp_path: Path, pools: Any
) -> None:
    job_id = make_job(sessions, tmp_path, 2)
    JobRunner(sessions, workers=1, executor_factory=pools).submit(job_id)
    first, future = pools.submitted[0]
    future.set_exception(BrokenProcessPool("a worker died"))
    eventually(lambda: len(pools.made) == 2 and len(pools.made[1].submitted) == 1)
    assert pools.made[0].shut_down
    assert len(pools.made[1].submitted) == 1
    assert statuses(sessions, job_id) == (JobStatus.RUNNING, [X, R])
    with sessions() as session:
        sample = session.get(Sample, first.sample_id)
        assert sample is not None
        assert sample.error == "worker process crashed"


def test_samples_left_running_are_run_again_on_start(
    sessions: sessionmaker[Session], tmp_path: Path, pools: Any
) -> None:
    job_id = make_job(sessions, tmp_path, 2, status=R)
    JobRunner(sessions, workers=2, executor_factory=pools).start()
    assert [task.name for task, _ in pools.submitted] == ["s0", "s1"]
    assert statuses(sessions, job_id) == (JobStatus.RUNNING, [R, R])


def test_shutdown_drops_queued_samples_and_records_nothing_more(
    sessions: sessionmaker[Session], tmp_path: Path, pools: Any
) -> None:
    job_id = make_job(sessions, tmp_path, 2)
    runner = JobRunner(sessions, workers=1, executor_factory=pools)
    runner.submit(job_id)
    runner.shutdown()
    task, future = pools.submitted[0]
    future.set_result(result_for(task))
    time.sleep(0.2)  # time for a recording thread to run, if one wrongly did
    assert len(pools.submitted) == 1
    assert pools.made[0].shut_down
    assert statuses(sessions, job_id) == (JobStatus.RUNNING, [R, Q])


def test_cancelling_drops_samples_not_started_and_lets_running_ones_finish(
    sessions: sessionmaker[Session], tmp_path: Path, pools: Any
) -> None:
    job_id = make_job(sessions, tmp_path, 3)
    runner = JobRunner(sessions, workers=1, executor_factory=pools)
    runner.submit(job_id)
    runner.cancel(job_id)
    assert statuses(sessions, job_id) == (JobStatus.CANCELLING, [R, C, C])

    task, future = pools.submitted[0]
    future.set_result(result_for(task))
    eventually(lambda: statuses(sessions, job_id)[0] == JobStatus.CANCELLED)
    # The running sample's result is kept, no other samples made it to the cpu pool
    assert statuses(sessions, job_id) == (JobStatus.CANCELLED, [F, C, C])
    assert len(pools.submitted) == 1
    with sessions() as session:
        job = session.get(Job, job_id)
        assert job is not None
        assert job.finished_at is not None


def test_cancelling_a_job_before_it_starts_stops_it_at_once(
    sessions: sessionmaker[Session], tmp_path: Path, pools: Any
) -> None:
    job_id = make_job(sessions, tmp_path, 2)
    runner = JobRunner(sessions, workers=1, executor_factory=pools)
    runner.cancel(job_id)
    assert statuses(sessions, job_id) == (JobStatus.CANCELLED, [C, C])
    runner.submit(job_id)
    assert pools.submitted == []


def test_a_restart_finishes_cancelling_rather_than_running_again(
    sessions: sessionmaker[Session], tmp_path: Path, pools: Any
) -> None:
    job_id = make_job(sessions, tmp_path, 2, status=R)
    with sessions() as session:
        job = session.get(Job, job_id)
        assert job is not None
        job.status = JobStatus.CANCELLING
        session.commit()
    JobRunner(sessions, workers=2, executor_factory=pools).start()
    assert pools.submitted == []
    assert statuses(sessions, job_id) == (JobStatus.CANCELLED, [C, C])


# A real process pool whose every task kills its worker process. It runs in a
# subprocess, so if the runner deadlocks the test fails instead of hanging pytest.
_CRASHING_POOL_RUN = textwrap.dedent(
    """
    import json, multiprocessing, os, sys, time
    from concurrent.futures import ProcessPoolExecutor
    from scalehd_server.db import make_engine
    from scalehd_server.models import Base, Job, Sample, User
    from scalehd_server.runner import JobRunner
    from scalehd_server.schemas import GenotypeMethod, JobSettings
    from sqlalchemy.orm import sessionmaker

    class CrashingPool(ProcessPoolExecutor):
        def submit(self, fn, /, *args, **kwargs):
            return super().submit(os._exit, 1)

    def main(folder):
        engine = make_engine(f"sqlite:///{folder}/crash.db")
        Base.metadata.create_all(engine)
        sessions = sessionmaker(engine)
        with sessions() as session:
            job = Job(
                owner=User(username="autotest-user", password_hash="x"),
                name="Demo",
                settings=JobSettings(method=GenotypeMethod.MODEL).model_dump(mode="json"),
                output_dir=f"{folder}/ws",
            )
            job.samples = [Sample(name=f"s{i}") for i in range(2)]
            session.add(job)
            session.commit()
            job_id = job.id
        context = multiprocessing.get_context("forkserver")
        runner = JobRunner(sessions, 1, lambda n: CrashingPool(n, mp_context=context))
        runner.submit(job_id)
        deadline = time.monotonic() + 20
        while time.monotonic() < deadline:
            with sessions() as session:
                if session.get(Job, job_id).status == "finished":
                    break
            time.sleep(0.1)
        with sessions() as session:
            job = session.get(Job, job_id)
            print(json.dumps([job.status, [[s.status, s.error] for s in job.samples]]), flush=True)
        os._exit(0)  # don't wait for a stuck pool thread at exit

    if __name__ == "__main__":
        main(sys.argv[1])
    """
)


def test_a_real_pool_whose_workers_crash_fails_each_sample_without_hanging(
    tmp_path: Path,
) -> None:
    script = tmp_path / "crash.py"
    script.write_text(_CRASHING_POOL_RUN)
    done = subprocess.run(
        [sys.executable, str(script), str(tmp_path)], capture_output=True, text=True, timeout=90
    )
    assert done.returncode == 0, done.stderr
    status, samples = json.loads(done.stdout.strip().splitlines()[-1])
    assert samples == [["failed", "worker process crashed"]] * 2
    assert status == "finished"
