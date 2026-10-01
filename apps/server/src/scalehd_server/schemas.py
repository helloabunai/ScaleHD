"""Request and response bodies of the HTTP API.

The web frontend mirrors these in ``apps/web/src/api.ts``.
"""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime
from typing import Annotated, Any, Literal

from pydantic import AfterValidator, BaseModel, ConfigDict, Field, StringConstraints
from scalehd.genotype import CallerSettings
from scalehd.pairs import DiscordancePolicy

from .models import JobStatus, SampleStatus


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


class JobSettings(BaseModel):
    """What a job runs, and with which settings. Unset thresholds keep the core defaults."""

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
    settings: JobSettings = JobSettings()


class SampleOut(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: int
    name: str
    r1: str
    r2: str | None
    status: SampleStatus
    genotype: str | None
    quality: float | None
    flags: list[str]
    error: str | None


class JobSummary(BaseModel):
    id: int
    name: str
    status: JobStatus
    created_at: datetime
    started_at: datetime | None
    finished_at: datetime | None
    sample_count: int
    samples_done: int


class JobOut(JobSummary):
    settings: JobSettings
    samples: list[SampleOut]
    error: str | None
