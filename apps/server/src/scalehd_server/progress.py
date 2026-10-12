"""Sample progress monitor. Dumps to json, deleted upon completion.
"""

from __future__ import annotations

import contextlib
import os
import time
from enum import StrEnum
from pathlib import Path

from pydantic import BaseModel, Field, ValidationError

FILE = "progress.json"


class SampleStage(StrEnum):

    SIMULATING = "simulating" # not for real jobs
    COUNTING = "counting"
    CALLING = "calling"


class SampleProgress(BaseModel):
    stage: SampleStage
    share: float | None = Field(default=None, ge=0, le=1)


def read_progress(folder: Path) -> SampleProgress | None:
    try:
        return SampleProgress.model_validate_json((folder / FILE).read_bytes())
    except (OSError, ValidationError):
        return None


def clear_progress(folder: Path) -> None:
    with contextlib.suppress(OSError):
        (folder / FILE).unlink(missing_ok=True)


class ProgressReporter:
    """writes progress/curr stage."""

    def __init__(self, folder: Path, every: float = 0.5) -> None:
        self._folder = folder
        self._every = every
        self._stage: SampleStage | None = None
        self._written = 0.0

    def stage(self, stage: SampleStage, share: float | None = None) -> None:
        now = time.monotonic()
        if stage == self._stage and now - self._written < self._every:
            return
        self._stage, self._written = stage, now
        self._write(SampleProgress(stage=stage, share=share))

    def counted(self, share: float) -> None:
        # For scalehd.count_fastq
        self.stage(SampleStage.COUNTING, share)

    def _write(self, progress: SampleProgress) -> None:
        self._folder.mkdir(parents=True, exist_ok=True)
        aside = self._folder / f".{FILE}.part"
        aside.write_text(progress.model_dump_json())
        os.replace(aside, self._folder / FILE)
