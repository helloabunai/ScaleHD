"""Each user's default settings for new jobs. Not built yet."""

from fastapi import APIRouter

from ..auth import CurrentUser
from ..db import DbSession
from ..errors import not_implemented
from ..schemas import JobSettings

router = APIRouter(prefix="/settings", tags=["settings"])


@router.get("")
def get_settings(user: CurrentUser) -> JobSettings:
    raise not_implemented("settings")


@router.put("")
def save_settings(new: JobSettings, user: CurrentUser, session: DbSession) -> JobSettings:
    raise not_implemented("settings")
