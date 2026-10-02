"""Each user's default settings for new jobs."""

from fastapi import APIRouter

from ..auth import CurrentUser
from ..db import DbSession
from ..schemas import JobSettings

router = APIRouter(prefix="/settings", tags=["settings"])


@router.get("")
def get_settings(user: CurrentUser) -> JobSettings:
    return JobSettings.model_validate(user.default_settings)


@router.put("")
def save_settings(new: JobSettings, user: CurrentUser, session: DbSession) -> JobSettings:
    user.default_settings = new.model_dump(mode="json")
    session.commit()
    return new
