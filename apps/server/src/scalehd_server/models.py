"""Database tables: accounts, login sessions, jobs, and the samples in each job."""

from __future__ import annotations

from datetime import UTC, datetime
from enum import StrEnum
from typing import Any

from sqlalchemy import JSON, ForeignKey, String
from sqlalchemy.orm import Mapped, mapped_column, relationship

from .db import Base


def _now() -> datetime:
    return datetime.now(UTC)


class JobStatus(StrEnum):
    QUEUED = "queued"
    RUNNING = "running"
    FINISHED = "finished"  # every sample has run, though some may have failed
    FAILED = "failed"  # the job itself could not run
    CANCELLED = "cancelled"


class Theme(StrEnum):
    """Light or dark pages, or whichever the user's computer is set to."""

    SYSTEM = "system"
    LIGHT = "light"
    DARK = "dark"


class SampleStatus(StrEnum):
    QUEUED = "queued"
    RUNNING = "running"
    FINISHED = "finished"
    FAILED = "failed"


class User(Base):
    __tablename__ = "users"

    id: Mapped[int] = mapped_column(primary_key=True)
    username: Mapped[str] = mapped_column(String(64), unique=True)
    password_hash: Mapped[str]
    is_admin: Mapped[bool] = mapped_column(default=False)
    created_at: Mapped[datetime] = mapped_column(default=_now)
    # Starting settings for this user's new jobs, as a schemas.JobSettings dict.
    default_settings: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    theme: Mapped[Theme] = mapped_column(default=Theme.SYSTEM)

    jobs: Mapped[list[Job]] = relationship(back_populates="owner")
    sessions: Mapped[list[LoginSession]] = relationship(
        back_populates="user", cascade="all, delete-orphan"
    )


class LoginSession(Base):
    """One logged-in browser. The cookie holds a random token, only its hash is kept."""

    __tablename__ = "login_sessions"

    id: Mapped[int] = mapped_column(primary_key=True)
    user_id: Mapped[int] = mapped_column(ForeignKey("users.id", ondelete="CASCADE"))
    token_hash: Mapped[str] = mapped_column(String(64), unique=True)
    created_at: Mapped[datetime] = mapped_column(default=_now)
    expires_at: Mapped[datetime]

    user: Mapped[User] = relationship(back_populates="sessions")


class Job(Base):
    """A named batch of samples run with one set of settings."""

    __tablename__ = "jobs"

    id: Mapped[int] = mapped_column(primary_key=True)
    owner_id: Mapped[int] = mapped_column(ForeignKey("users.id"))
    name: Mapped[str] = mapped_column(String(200))
    status: Mapped[JobStatus] = mapped_column(default=JobStatus.QUEUED)
    settings: Mapped[dict[str, Any]] = mapped_column(JSON)
    created_at: Mapped[datetime] = mapped_column(default=_now)
    started_at: Mapped[datetime | None]
    finished_at: Mapped[datetime | None]
    error: Mapped[str | None]
    demo: Mapped[bool] = mapped_column(default=False)
    # The job's folder in the workspace.
    output_dir: Mapped[str | None]

    owner: Mapped[User] = relationship(back_populates="jobs")
    samples: Mapped[list[Sample]] = relationship(
        back_populates="job", cascade="all, delete-orphan", order_by="Sample.id"
    )


class Sample(Base):
    """One sample's input and, once run, its genotype call.

    Counts and call files are in the sample's folder inside the job's folder. The call
    is also stored here, with its label, confidence and flags copied out so job pages can
    list them without reading the JSON.
    """

    __tablename__ = "samples"

    id: Mapped[int] = mapped_column(primary_key=True)
    job_id: Mapped[int] = mapped_column(ForeignKey("jobs.id"))
    name: Mapped[str] = mapped_column(String(200))
    # Paths of real input files, inside the data root. None for simulated samples.
    r1: Mapped[str | None]
    r2: Mapped[str | None]
    # Simulated samples: {"alleles": [labels], "pairs": n, "seed": s}.
    simulation: Mapped[dict[str, Any] | None] = mapped_column(JSON)
    # Simulated samples: the true genotype label, and whether the call matched it.
    truth: Mapped[str | None]
    matches_truth: Mapped[bool | None]
    status: Mapped[SampleStatus] = mapped_column(default=SampleStatus.QUEUED)
    genotype: Mapped[str | None]
    confidence: Mapped[float | None]
    flags: Mapped[list[str]] = mapped_column(JSON, default=list)
    call: Mapped[dict[str, Any] | None] = mapped_column(JSON)  # scalehd.call/1
    error: Mapped[str | None]

    job: Mapped[Job] = relationship(back_populates="samples")
