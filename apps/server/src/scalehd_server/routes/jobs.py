"""List and follow jobs. Spawn jobs/fetch each sample's results.

The web interface follows progress by polling ``GET /jobs/{job_id}``. Server-sent
events could replace polling later without changing anything else.
"""

from pathlib import Path

from fastapi import APIRouter, HTTPException, Response, status
from fastapi.responses import FileResponse
from sqlalchemy import select
from sqlalchemy.orm import selectinload

from .. import demo
from ..auth import CurrentUser
from ..config import ServerConfig
from ..db import DbSession
from ..errors import not_implemented
from ..models import Job, JobStatus, Sample, User
from ..results import FileKind, sample_detail, sample_file
from ..runner import Runner
from ..schemas import (
    JobCreate,
    JobOut,
    JobSettings,
    JobSummary,
    SampleDetail,
    job_out,
    job_summary,
)
from ..worker import RUNNABLE_METHODS
from ..workspace import WorkspaceError, remove_job_folder, write_job_folder

router = APIRouter(prefix="/jobs", tags=["jobs"])


@router.get("")
def list_jobs(user: CurrentUser, session: DbSession) -> list[JobSummary]:
    """The user's own jobs, newest first."""
    jobs = session.scalars(
        select(Job)
        .where(Job.owner_id == user.id)
        .options(selectinload(Job.samples))
        .order_by(Job.created_at.desc(), Job.id.desc())
    ).all()
    return [job_summary(job) for job in jobs]


@router.post("", status_code=status.HTTP_201_CREATED)
def create_job(
    job: JobCreate, user: CurrentUser, session: DbSession, config: ServerConfig, runner: Runner
) -> JobOut:
    """Check each sample's files exist in the input directory, store the job, queue it."""
    settings = job.settings or JobSettings.model_validate(user.default_settings)
    if settings.method not in RUNNABLE_METHODS:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            f"{settings.method} genotyping is not available yet. Placeholder flag.",
        )
    raise not_implemented("jobs")


@router.post("/demo", status_code=status.HTTP_201_CREATED)
def create_demo_job(
    user: CurrentUser, session: DbSession, config: ServerConfig, runner: Runner
) -> JobOut:
    """Simulated samples with known genotypes, run like any other job."""
    job = demo.demo_job(user, JobSettings.model_validate(user.default_settings))
    session.add(job)
    session.flush()
    try:
        job.output_dir = str(write_job_folder(config.workspace, job))
    except WorkspaceError as exc:
        session.rollback()
        raise HTTPException(status.HTTP_500_INTERNAL_SERVER_ERROR, str(exc)) from None
    session.commit()
    runner.submit(job.id)
    return job_out(job)


def _own_job(session: DbSession, user: User, job_id: int) -> Job:
    job = session.get(Job, job_id)
    # Someone else's job looks the same as a missing one, so ids can't be probed.
    if job is None or job.owner_id != user.id:
        raise HTTPException(status.HTTP_404_NOT_FOUND, "no such job")
    return job


@router.get("/{job_id}")
def get_job(job_id: int, user: CurrentUser, session: DbSession) -> JobOut:
    return job_out(_own_job(session, user, job_id))


@router.post("/{job_id}/cancel")
def cancel_job(job_id: int, user: CurrentUser, session: DbSession, runner: Runner) -> JobOut:
    raise not_implemented("jobs")


@router.delete("/{job_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_job(job_id: int, user: CurrentUser, session: DbSession, config: ServerConfig) -> None:
    """Delete a job that has stopped (job folder in workspace, and DB records)
    A running job is refused. cos it's running.
    """
    job = _own_job(session, user, job_id)
    if job.status in (JobStatus.QUEUED, JobStatus.RUNNING):
        raise HTTPException(status.HTTP_409_CONFLICT, "the job is still running")
    if job.output_dir is not None:
        try:
            remove_job_folder(config.workspace, user.username, Path(job.output_dir))
        except WorkspaceError as exc:
            raise HTTPException(status.HTTP_500_INTERNAL_SERVER_ERROR, str(exc)) from None
    session.delete(job)
    session.commit()


def _job_sample(job: Job, sample_id: int) -> Sample:
    for sample in job.samples:
        if sample.id == sample_id:
            return sample
    raise HTTPException(status.HTTP_404_NOT_FOUND, "no such sample")


@router.get("/{job_id}/samples/{sample_id}")
def get_sample(job_id: int, sample_id: int, user: CurrentUser, session: DbSession) -> SampleDetail:
    """Sample results = genotype call, read distributions and summary info.
    Files for download not finalised."""
    job = _own_job(session, user, job_id)
    return sample_detail(job, _job_sample(job, sample_id))


@router.get("/{job_id}/samples/{sample_id}/files/{kind}", response_class=FileResponse)
def get_sample_file(
    job_id: int, sample_id: int, kind: FileKind, user: CurrentUser, session: DbSession
) -> FileResponse:
    """Download sample results + associated data"""
    job = _own_job(session, user, job_id)
    sample = _job_sample(job, sample_id)
    path = sample_file(job, sample, kind)
    if path is None or not path.is_file():
        raise HTTPException(status.HTTP_404_NOT_FOUND, "no such file")
    name = path.name if kind in ("r1", "r2") else f"{path.parent.name}-{path.name}"
    return FileResponse(path, filename=name)


@router.get("/{job_id}/report", response_class=Response)
def job_report(job_id: int, user: CurrentUser, session: DbSession) -> Response:
    """A PDF summary of the job's calls and flags."""
    raise not_implemented("reports")
