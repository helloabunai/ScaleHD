"""Hands samples to worker processes and records what comes back.

One first-come, first-served queue of samples across all jobs. At most ``workers``
samples are with the pool at once, and only those are marked running. Each result is
saved as its sample finishes. Only this process writes the database.
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
from typing import Annotated

from fastapi import Depends, Request
from sqlalchemy import select, update
from sqlalchemy.orm import Session, sessionmaker

from .models import Job, JobStatus, Sample, SampleStatus
from .schemas import JobSettings
from .worker import SampleResult, SampleTask, run_task
from .workspace import sample_folders

ExecutorFactory = Callable[[int], Executor]

_DONE = (SampleStatus.FINISHED, SampleStatus.FAILED)


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
    ) -> None:
        self._sessions = sessions
        self._workers = workers
        self._executor_factory = executor_factory
        self._pool: Executor | None = None
        self._queue: deque[int] = deque()
        # Pooled samples and if a lock is present.
        self._handed_out: dict[Future[SampleResult], tuple[int, Executor]] = {}
        self._lock = threading.RLock()
        self._stopped = False

    def start(self) -> None:
        """Pick up where the last run stopped. Samples left running are re-queued."""
        with self._sessions() as session:
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
            while not self._stopped and self._queue and len(self._handed_out) < self._workers:
                task = self._start_sample(self._queue.popleft())
                if task is None:
                    continue
                pool = self._executor()
                future = pool.submit(run_task, task)
                self._handed_out[future] = (task.sample_id, pool)
                future.add_done_callback(self._finished_soon)

    def _finished_soon(self, future: Future[SampleResult]) -> None:
        threading.Thread(target=self._finished, args=(future,), daemon=True).start()

    def _finished(self, future: Future[SampleResult]) -> None:
        with self._lock:
            entry = self._handed_out.pop(future, None)
            if entry is None or self._stopped:
                return
            sample_id, pool = entry
            result: SampleResult | None = None
            error: str | None = None
            try:
                result = future.result()
            except CancelledError:
                return
            except BrokenProcessPool:
                error = "worker process crashed"
                self._discard(pool)
            except Exception as exc:
                error = f"{type(exc).__name__}: {exc}"
            self._record(sample_id, result, error)
        self._fill()

    def _start_sample(self, sample_id: int) -> SampleTask | None:
        """Mark a queued sample running (and its job, if first) and describe it."""
        with self._sessions() as session:
            sample = session.get(Sample, sample_id)
            if sample is None or sample.status != SampleStatus.QUEUED:
                return None
            job = sample.job
            sample.status = SampleStatus.RUNNING
            if job.status == JobStatus.QUEUED:
                job.status, job.started_at = JobStatus.RUNNING, datetime.now(UTC)
            task = SampleTask(
                sample_id=sample.id,
                name=sample.name,
                folder=sample_folders(job)[sample.id],
                settings=JobSettings.model_validate(job.settings),
                r1=Path(sample.r1) if sample.r1 else None,
                r2=Path(sample.r2) if sample.r2 else None,
                simulation=sample.simulation,
                truth=sample.truth,
            )
            session.commit()
            return task

    def _record(self, sample_id: int, result: SampleResult | None, error: str | None) -> None:
        with self._sessions() as session:
            sample = session.get(Sample, sample_id)
            if sample is None:
                return
            if result is None:
                sample.status, sample.error = SampleStatus.FAILED, error
            else:
                sample.status = SampleStatus.FINISHED
                sample.call, sample.genotype = result.call, result.genotype
                sample.confidence, sample.flags = result.confidence, result.flags
                sample.matches_truth = result.matches_truth
            job = sample.job
            if all(other.status in _DONE for other in job.samples):
                job.status, job.finished_at = JobStatus.FINISHED, datetime.now(UTC)
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


def get_runner(request: Request) -> JobRunner:
    runner: JobRunner = request.app.state.runner
    return runner


Runner = Annotated[JobRunner, Depends(get_runner)]
