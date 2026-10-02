"""Where sequencing data comes from and where results go, so users can find them."""

from fastapi import APIRouter

from ..auth import CurrentUser
from ..config import ServerConfig
from ..schemas import Folders
from ..workspace import user_folder

router = APIRouter(prefix="/folders", tags=["folders"])


@router.get("")
def folders(user: CurrentUser, config: ServerConfig) -> Folders:
    return Folders(
        data_root=str(config.data_root) if config.data_root else None,
        workspace=str(config.workspace),
        your_folder=str(user_folder(config.workspace, user.username)),
    )
