"""Hands samples to worker processes and records what comes back.

One first-come, first-served queue of samples across all jobs. Samples share ``workers``
cores. A legacy job's sample takes its (user settings) cores per sample. A sample
waits at the front of the queue until enough are free. Only samples with the pool are
marked running. Each result is saved as its sample finishes. Only this process writes the
database.
"""

from __future__ import annotations

import multiprocessing
import threading
from collections import deque
from collections.abc import Callable
from concurrent.futures import CancelledError, Executor, Future, ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from datetime import UTC, datetime
from pathlib import Path
from typing import Annotated, Literal

from fastapi import Depends, Request
from sqlalchemy import select, update
from sqlalchemy.orm import Session, sessionmaker

from .legacy import LegacyError
from .models import Job, JobStatus, Sample, SampleStatus
from .progress import clear_progress
from .schemas import GenotypeMethod, JobSettings
from .worker import SampleResult, SampleTask, run_task
from .workspace import sample_folders

ExecutorFactory = Callable[[int], Executor]

# "final" statuses i.e. won't change from these once set
_DONE = (SampleStatus.FINISHED, SampleStatus.FAILED, SampleStatus.CANCELLED)


def process_pool(workers: int) -> Executor:
    # forkserver, because forking a process that already runs threads (uvicorn's) is
    # unsafe.
    return ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("forkserver"))


class JobRunner:
    def __init__(
        self,
        sessions: sessionmaker[Session],
        workers: int,
        executor_factory: ExecutorFactory = process_pool,
        legacy_runner: Path | None = None,
    ) -> None:
        self._sessions = sessions
        self._workers = workers
        self._executor_factory = executor_factory
        self._legacy_runner = legacy_runner
        self._pool: Executor | None = None
        self._queue: deque[int] = deque()
        # pooled sample id, pool and core(s).
        self._handed_out: dict[Future[SampleResult], tuple[int, Executor, int]] = {}
        self._lock = threading.RLock()
        self._stopped = False

    def start(self) -> None:
        """Pick up where the last run stopped. Samples left running are re-queued, unless
        their job was being cancelled."""
        with self._sessions() as session:
            cancelling = session.scalars(
                select(Job).where(Job.status == JobStatus.CANCELLING)
            ).all()
            for job in cancelling:
                for sample in job.samples:
                    if sample.status not in _DONE:
                        sample.status = SampleStatus.CANCELLED
                _update_job(job)
            session.execute(
                update(Sample)
                .where(Sample.status == SampleStatus.RUNNING)
                .values(status=SampleStatus.QUEUED)
            )
            session.commit()
            queued = session.scalars(
                select(Sample.id)
                .join(Job)
                .where(Sample.status == SampleStatus.QUEUED)
                .order_by(Job.created_at, Job.id, Sample.id)
            ).all()
        with self._lock:
            self._queue.extend(queued)
        self._fill()

    def submit(self, job_id: int) -> None:
        """Queue a stored job's queued samples."""
        with self._sessions() as session:
            queued = session.scalars(
                select(Sample.id)
                .where(Sample.job_id == job_id, Sample.status == SampleStatus.QUEUED)
                .order_by(Sample.id)
            ).all()
        with self._lock:
            self._queue.extend(queued)
        self._fill()

    def cancel(self, job_id: int) -> None:
        """Stop a queued or running job. Job samples not yet started won't be.
        Those already in the processor pool finish and are recorded as usual, and the job is
        cancelled once they have (messing with pools is bad smell)."""
        with self._lock:
            with self._sessions() as session:
                job = session.get(Job, job_id)
                if job is None or job.status not in (JobStatus.QUEUED, JobStatus.RUNNING):
                    return
                dropped = set()
                for sample in job.samples:
                    if sample.status == SampleStatus.QUEUED:
                        sample.status = SampleStatus.CANCELLED
                        dropped.add(sample.id)
                job.status = JobStatus.CANCELLING
                _update_job(job)
                session.commit()
            self._queue = deque(i for i in self._queue if i not in dropped)

    def shutdown(self) -> None:
        """Stop handing out samples and drop those not started.

        Samples already with the pool stay marked running, so the next start queues
        them again. Nothing more is recorded.
        """
        with self._lock:
            self._stopped = True
            self._queue.clear()
            pool, self._pool = self._pool, None
        if pool is not None:
            pool.shutdown(wait=False, cancel_futures=True)

    def _fill(self) -> None:
        with self._lock:
            while not self._stopped and self._queue:
                busy = sum(cores for _, _, cores in self._handed_out.values())
                task = self._start_sample(self._queue[0], self._workers - busy)
                if task == "wait":
                    break
                self._queue.popleft()
                if task is None:
                    continue
                pool = self._executor()
                future = pool.submit(run_task, task)
                self._handed_out[future] = (task.sample_id, pool, task.threads)
                future.add_done_callback(self._finished_soon)

    def _finished_soon(self, future: Future[SampleResult]) -> None:
        threading.Thread(target=self._finished, args=(future,), daemon=True).start()

    def _finished(self, future: Future[SampleResult]) -> None:
        with self._lock:
            entry = self._handed_out.pop(future, None)
            if entry is None or self._stopped:
                return
            sample_id, pool, _ = entry
            result: SampleResult | None = None
            error: str | None = None
            try:
                result = future.result()
            except CancelledError:
                return
            except BrokenProcessPool:
                error = "worker process crashed"
                self._discard(pool)
            except LegacyError as exc:
                # Already says it was ScaleHD 1.x, and where it failed.
                error = str(exc)
            except Exception as exc:
                error = f"{type(exc).__name__}: {exc}"
            self._record(sample_id, result, error)
        self._fill()

    def _start_sample(self, sample_id: int, free: int) -> SampleTask | Literal["wait"] | None:
        """Mark a queued sample running (and its job, if first) and describe it. "wait" if
        it needs more cores than are free, None if it isn't queued any more."""
        with self._sessions() as session:
            sample = session.get(Sample, sample_id)
            if sample is None or sample.status != SampleStatus.QUEUED:
                return None
            job = sample.job
            settings = JobSettings.model_validate(job.settings)
            threads = self._cores(settings)
            if threads > free:
                return "wait"
            sample.status = SampleStatus.RUNNING
            if job.status == JobStatus.QUEUED:
                job.status, job.started_at = JobStatus.RUNNING, datetime.now(UTC)
            folder = sample_folders(job)[sample.id]
            # Left by an earlier attempt/cancelled job/easier to clean + go fresh
            clear_progress(folder)
            task = SampleTask(
                sample_id=sample.id,
                name=sample.name,
                folder=folder,
                settings=settings,
                r1=Path(sample.r1) if sample.r1 else None,
                r2=Path(sample.r2) if sample.r2 else None,
                simulation=sample.simulation,
                truth=sample.truth,
                threads=threads,
                legacy_runner=self._legacy_runner,
            )
            session.commit()
            return task

    def _cores(self, settings: JobSettings) -> int:
        """Cores for one of a job's samples"""
        if settings.method is GenotypeMethod.LEGACY and settings.legacy is not None:
            return min(settings.legacy.cores_per_sample, self._workers)
        return 1

    def _record(self, sample_id: int, result: SampleResult | None, error: str | None) -> None:
        with self._sessions() as session:
            sample = session.get(Sample, sample_id)
            if sample is None:
                return
            clear_progress(sample_folders(sample.job)[sample.id])
            if result is None:
                sample.status, sample.error = SampleStatus.FAILED, error
            else:
                sample.status = SampleStatus.FINISHED
                sample.call, sample.genotype = result.call, result.genotype
                sample.confidence, sample.flags = result.confidence, result.flags
                sample.matches_truth = result.matches_truth
            _update_job(sample.job)
            session.commit()

    def _executor(self) -> Executor:
        # Made on first use, so the server starts instantly.
        if self._pool is None:
            self._pool = self._executor_factory(self._workers)
        return self._pool

    def _discard(self, pool: Executor) -> None:
        # A crashed worker breaks its whole pool. Every sampled queued
        # it held fails with BrokenProcessPool
        if self._pool is pool:
            self._pool = None
            pool.shutdown(wait=False)


def _update_job(job: Job) -> None:
    """Finish a job whose samples are all done."""
    if any(sample.status not in _DONE for sample in job.samples):
        return
    cancelled = job.status == JobStatus.CANCELLING
    job.status = JobStatus.CANCELLED if cancelled else JobStatus.FINISHED
    job.finished_at = datetime.now(UTC)


def get_runner(request: Request) -> JobRunner:
    runner: JobRunner = request.app.state.runner
    return runner


Runner = Annotated[JobRunner, Depends(get_runner)]
