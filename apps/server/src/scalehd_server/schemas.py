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

from .models import Job, JobStatus, SampleStatus, Theme


class Health(BaseModel):
    status: Literal["ok"]
    version: str
    core_version: str
    python: str
    platform: str
    sqlite: str
    libraries: dict[str, str]


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
    theme: Theme


class ThemeChange(BaseModel):
    theme: Theme


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

    # Model-based default until the legacy method exists as an actual optino
    method: GenotypeMethod = GenotypeMethod.MODEL
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
    """A sample's FASTQ files, as paths relative to the server's data rpot."""

    name: str
    r1: str
    r2: str | None = None


class InputSample(BaseModel):
    """A sample found in the data folder: its FASTQ files, paired by name."""

    name: str
    files: list[str]
    r1: str | None  # None when it can't be run for any reason
    r2: str | None  # None when it can't be run for any reason
    size: int
    undetermined: bool  # couldn't regex the sample name from filename.
    skipped: str | None  # reason for skip e.g. R2 without R1


class InputSubfolder(BaseModel):
    """A folder inside the open one, with what it holds, for the folder tree."""

    name: str
    folders: int
    samples: int  # runnable samples (name regex OK)


class InputFolder(BaseModel):
    """a run/data dir in the server's data root."""

    folder: str  # Relative to the data folder; "" is the data folder itself.
    path: str  # complete dir path
    folders: list[InputSubfolder]
    samples: list[InputSample]
    other_files: list[str]  # non FASTQ files if any


MAX_TAGS = 5
MAX_TAG_LENGTH = 15
TagName = Annotated[
    str, StringConstraints(strip_whitespace=True, min_length=1, max_length=MAX_TAG_LENGTH)
]


class JobTag(BaseModel):
    """self explanatory"""

    model_config = ConfigDict(from_attributes=True)

    id: int
    name: str


class TagOut(JobTag):
    """A tag, with how many jobs use it (from any user)."""

    jobs: int


class TagChange(BaseModel):
    """Making or renaming a tag."""

    name: TagName


class JobTags(BaseModel):
    """Tags sorted by id."""

    tags: list[int] = Field(max_length=MAX_TAGS)


class JobCreate(BaseModel):
    name: str = Field(min_length=1, max_length=200)
    samples: list[InputPair] = Field(min_length=1)
    settings: JobSettings | None = None  # None = unset = default
    tags: list[int] = Field(default_factory=list, max_length=MAX_TAGS)


class AdminChange(BaseModel):
    is_admin: bool


class SampleOut(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: int
    name: str
    r1: str | None
    r2: str | None
    status: SampleStatus
    genotype: str | None
    confidence: float | None
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
    tags: list[JobTag]


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
        "tags": [JobTag.model_validate(tag) for tag in job.tags],
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


class CagBar(BaseModel):
    cag: int
    # Molecules at this length, including those read to the tract's end too near the
    # read's own end to confirm it.
    molecules: int
    # Molecules whose CAG is only known to be at least this (reads ended in the tract).
    lower_bound: int


class CagChart(BaseModel):
    """Molecules at each CAG length within one structure of the called alleles."""

    caacag: int
    ccgcca: int
    ccg: int
    cct: int
    # The called alleles on this structure, and their CAG lengths (to highlight).
    alleles: list[str]
    called: list[int]
    bars: list[CagBar]


class CcgBar(BaseModel):
    ccg: int
    molecules: int


class Cell(BaseModel):
    cag: int
    ccg: int
    molecules: int


class Reads(BaseModel):
    """Sample read information"""

    molecules: int
    complete: int
    partial: int
    dropped: int
    unusable: int
    read_outcomes: dict[str, int]
    discordant: dict[str, int]


class SampleDetail(BaseModel):
    """All data needed in an indv sample result page"""

    sample: SampleOut
    job_id: int
    job_name: str
    demo: bool
    tags: list[JobTag]
    folder: str | None
    # The full call (scalehd.call/2), once the sample has been called.
    call: dict[str, Any] | None
    cag_charts: list[CagChart]
    ccg: list[CcgBar]
    cells: list[Cell]
    reads: Reads | None
    # Which files can be downloaded: "call", "counts", "r1", "r2".
    files: list[str]
    previous_id: int | None
    next_id: int | None
