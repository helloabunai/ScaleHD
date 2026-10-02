"""Request and response bodies of the HTTP API.

The web frontend mirrors these in ``apps/web/src/api.ts``.
"""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime
from enum import StrEnum
from typing import Annotated, Any, Literal

from pydantic import AfterValidator, BaseModel, ConfigDict, Field, StringConstraints
from scalehd.genotype import CallerSettings
from scalehd.pairs import DiscordancePolicy

from .models import Job, JobStatus, SampleStatus


class Health(BaseModel):
    status: Literal["ok"]
    version: str
    core_version: str


# Usernames are case-insensitive, so they are stored and compared in lower case.
# Scope of users is like 4 people in a lab so these are not as strict requirements as
# a "real" webserver
_Username = Annotated[
    str, StringConstraints(strip_whitespace=True, max_length=64), AfterValidator(str.lower)
]
_NewPassword = Annotated[str, Field(min_length=8, max_length=256)]


class Login(BaseModel):
    username: _Username
    password: str = Field(max_length=256)


class NewAccount(BaseModel):
    # Letters, digits, dots, dashes and underscores, starting with a letter or digit.
    username: Annotated[_Username, StringConstraints(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]*$")]
    password: _NewPassword


class PasswordChange(BaseModel):
    current_password: str = Field(max_length=256)
    new_password: _NewPassword


class Folders(BaseModel):
    """Data folders from host machine i.e. sequencing data to pick from, and where results go."""

    data_root: str | None
    workspace: str
    # <workspace>/<username>: this user's jobs are saved here.
    your_folder: str


class Registration(BaseModel):
    """Whether accounts can be created, and whether the next one is the first (the admin)."""

    open: bool
    first_account: bool


class UserOut(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: int
    username: str
    is_admin: bool
    created_at: datetime


class GenotypeMethod(StrEnum):
    """How a job genotypes its samples. The web interface holds the names people see."""

    # ScaleHD 1.x: align reads to the reference library, then the 1.x genotyper.
    # Not runnable yet as it needs extracting from legacy/, alignment included.
    LEGACY = "legacy"
    # Read the repeat structure straight from each read, then the model-based caller
    # (scalehd.genotype). Work in progress.
    MODEL = "model"


class JobSettings(BaseModel):
    """What a job runs, and with which settings. Unset thresholds keep the core defaults."""

    method: GenotypeMethod = GenotypeMethod.LEGACY
    # Call genotypes, or only count each sample's repeat structures.
    call: bool = True
    discordant: DiscordancePolicy = DiscordancePolicy.DROP
    # Flag thresholds, as in scalehd.genotype.CallerSettings.
    min_molecules: int | None = Field(default=None, ge=0)
    min_posterior: float | None = Field(default=None, ge=0, le=1)
    max_background: float | None = Field(default=None, ge=0, le=1)
    max_dropped: float | None = Field(default=None, ge=0, le=1)

    def caller_settings(self) -> CallerSettings:
        thresholds = ("min_molecules", "min_posterior", "max_background", "max_dropped")
        overrides: dict[str, Any] = {
            name: value for name in thresholds if (value := getattr(self, name)) is not None
        }
        return replace(CallerSettings(), **overrides)


class InputPair(BaseModel):
    """A sample's FASTQ files, as paths relative to the server's input directory."""

    name: str
    r1: str
    r2: str | None = None


class JobCreate(BaseModel):
    name: str = Field(min_length=1, max_length=200)
    samples: list[InputPair] = Field(min_length=1)
    # Unset means the user's default settings.
    settings: JobSettings | None = None


class SampleOut(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: int
    name: str
    r1: str | None
    r2: str | None
    status: SampleStatus
    genotype: str | None
    quality: float | None
    flags: list[str]
    truth: str | None
    matches_truth: bool | None
    error: str | None


class JobSummary(BaseModel):
    id: int
    name: str
    demo: bool
    method: GenotypeMethod
    status: JobStatus
    created_at: datetime
    started_at: datetime | None
    finished_at: datetime | None
    output_dir: str | None
    sample_count: int
    samples_done: int


class JobOut(JobSummary):
    settings: JobSettings
    samples: list[SampleOut]
    error: str | None


def _summary(job: Job) -> dict[str, Any]:
    done = (SampleStatus.FINISHED, SampleStatus.FAILED)
    return {
        "id": job.id,
        "name": job.name,
        "demo": job.demo,
        "method": JobSettings.model_validate(job.settings).method,
        "status": job.status,
        "created_at": job.created_at,
        "started_at": job.started_at,
        "finished_at": job.finished_at,
        "output_dir": job.output_dir,
        "sample_count": len(job.samples),
        "samples_done": sum(sample.status in done for sample in job.samples),
    }


def job_summary(job: Job) -> JobSummary:
    """A job for the jobs list. Call inside an open session."""
    return JobSummary(**_summary(job))


def job_out(job: Job) -> JobOut:
    """A job with its samples. Call inside an open session."""
    return JobOut(
        **_summary(job),
        settings=JobSettings.model_validate(job.settings),
        samples=[SampleOut.model_validate(sample) for sample in job.samples],
        error=job.error,
    )
