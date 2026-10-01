"""Jobs: create, list, follow and cancel, and fetch each sample's results. Not built yet.

The web interface follows progress by polling ``GET /jobs/{job_id}``. Server-sent
events could replace polling later without changing anything else.
"""

from typing import Any

from fastapi import APIRouter, Response, status
from fastapi.responses import FileResponse

from ..auth import CurrentUser
from ..config import ServerConfig
from ..db import DbSession
from ..errors import not_implemented
from ..runner import Runner
from ..schemas import JobCreate, JobOut, JobSummary

router = APIRouter(prefix="/jobs", tags=["jobs"])


@router.get("")
def list_jobs(user: CurrentUser, session: DbSession) -> list[JobSummary]:
    """The user's jobs, newest first. Admins see everyone's."""
    raise not_implemented("jobs")


@router.post("", status_code=status.HTTP_201_CREATED)
def create_job(
    job: JobCreate, user: CurrentUser, session: DbSession, config: ServerConfig, runner: Runner
) -> JobOut:
    """Check each sample's files exist in the input directory, store the job, queue it."""
    raise not_implemented("jobs")


@router.get("/{job_id}")
def get_job(job_id: int, user: CurrentUser, session: DbSession) -> JobOut:
    raise not_implemented("jobs")


@router.post("/{job_id}/cancel")
def cancel_job(job_id: int, user: CurrentUser, session: DbSession, runner: Runner) -> JobOut:
    raise not_implemented("jobs")


@router.delete("/{job_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_job(job_id: int, user: CurrentUser, session: DbSession) -> None:
    """Delete a finished or cancelled job and its output files."""
    raise not_implemented("jobs")


@router.get("/{job_id}/samples/{sample_id}/call")
def sample_call(
    job_id: int, sample_id: int, user: CurrentUser, session: DbSession
) -> dict[str, Any]:
    """The full genotype call (``scalehd.call/1`` JSON)."""
    raise not_implemented("jobs")


@router.get("/{job_id}/samples/{sample_id}/counts", response_class=FileResponse)
def sample_counts(
    job_id: int, sample_id: int, user: CurrentUser, config: ServerConfig, session: DbSession
) -> FileResponse:
    """The sample's molecule counts (``scalehd.counts/1`` JSON), as written by the runner."""
    raise not_implemented("jobs")


@router.get("/{job_id}/report", response_class=Response)
def job_report(job_id: int, user: CurrentUser, session: DbSession) -> Response:
    """A PDF summary of the job's calls and flags."""
    raise not_implemented("reports")
