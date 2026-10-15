"""Each user's default settings for new jobs, and what this server can run jobs with."""

from fastapi import APIRouter

from ..auth import CurrentUser
from ..config import ServerConfig
from ..db import DbSession
from ..schemas import JobSettings, ServerOptions
from ..worker import available_methods

router = APIRouter(prefix="/settings", tags=["settings"])


@router.get("")
def get_settings(user: CurrentUser) -> JobSettings:
    """With ScaleHD 1.x's settings filled in, so they can be edited before they're used."""
    return JobSettings.model_validate(user.default_settings).with_legacy()


@router.put("")
def save_settings(new: JobSettings, user: CurrentUser, session: DbSession) -> JobSettings:
    user.default_settings = new.model_dump(mode="json")
    session.commit()
    return new.with_legacy()


@router.get("/options")
def get_options(user: CurrentUser, config: ServerConfig) -> ServerOptions:
    return ServerOptions(methods=available_methods(config), cores=config.workers)
