"""Runs jobs: each sample is counted and called in a worker process.

The runner lives in the server process. Workers are given file paths and settings
and return results, and only the server process writes to the database, so SQLite's
single writer is never contended. ``run_sample`` is the real work; queueing,
progress and cancelling are still stubs.
"""

from __future__ import annotations

import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Annotated, Any

from fastapi import Depends, Request
from scalehd.counts import count_fastq
from scalehd.genotype import call_genotype
from sqlalchemy.orm import Session, sessionmaker

from .schemas import GenotypeMethod, JobSettings

# Methods jobs can use. SHD 1.x joins once it is extracted from legacy/.
RUNNABLE_METHODS = frozenset({GenotypeMethod.MODEL})


def run_sample(
    r1: Path, r2: Path | None, out_dir: Path, settings: JobSettings
) -> dict[str, Any] | None:
    """Count one sample's molecules and call its genotype. Runs in a worker process.

    Writes ``counts.json`` and, when calling, ``call.json`` to out_dir. Returns the
    call, or None for a count-only job.
    """
    if settings.method not in RUNNABLE_METHODS:
        raise ValueError(f"{settings.method} genotyping is not available yet")
    out_dir.mkdir(parents=True, exist_ok=True)
    counts = count_fastq(r1, r2, policy=settings.discordant)
    counts.write_json(out_dir / "counts.json")
    if not settings.call:
        return None
    call = call_genotype(counts, settings.caller_settings()).to_dict()
    (out_dir / "call.json").write_text(json.dumps(call, indent=2) + "\n")
    return call


class JobRunner:
    def __init__(self, sessions: sessionmaker[Session], output_dir: Path, workers: int) -> None:
        self._sessions = sessions
        self._output_dir = output_dir
        self._workers = workers
        self._pool: ProcessPoolExecutor | None = None

    def start(self) -> None:
        """Pick up where the last run stopped.

        TODO: set samples left RUNNING by a crash or restart back to QUEUED, then
        submit every QUEUED job.
        """

    def submit(self, job_id: int) -> None:
        """Queue a stored job's samples.

        TODO: one ``run_sample`` per sample on the pool. As each finishes, record its
        status, genotype, quality, flags and call (or error); when the last one does,
        mark the job FINISHED.
        """
        raise NotImplementedError

    def cancel(self, job_id: int) -> None:
        """TODO: drop the job's queued samples and mark it CANCELLED once running ones end."""
        raise NotImplementedError

    def shutdown(self) -> None:
        if self._pool is not None:
            self._pool.shutdown(cancel_futures=True)
            self._pool = None

    def _executor(self) -> ProcessPoolExecutor:
        # Created on first use so the server starts instantly. forkserver, because
        # forking a process that already runs threads (uvicorn's) is unsafe.
        if self._pool is None:
            context = multiprocessing.get_context("forkserver")
            self._pool = ProcessPoolExecutor(self._workers, mp_context=context)
        return self._pool


def get_runner(request: Request) -> JobRunner:
    runner: JobRunner = request.app.state.runner
    return runner


Runner = Annotated[JobRunner, Depends(get_runner)]
