"""Browse the server's data folder for FASTQ files to run.

Files are under ``SCALEHD_DATA_ROOT``, mounted read-only from host to docker.
"""

from fastapi import APIRouter, HTTPException, status

from ..auth import CurrentUser
from ..config import ServerConfig
from ..inputs import InputError, list_folder
from ..schemas import InputFolder

router = APIRouter(prefix="/inputs", tags=["inputs"])


@router.get("")
def list_inputs(user: CurrentUser, config: ServerConfig, folder: str = "") -> InputFolder:
    """Subdir of the server's data root."""
    if config.data_root is None:
        raise HTTPException(
            status.HTTP_409_CONFLICT, "this server has no data folder set (SCALEHD_DATA_ROOT)"
        )
    try:
        return list_folder(config.data_root, folder)
    except InputError as exc:
        raise HTTPException(status.HTTP_404_NOT_FOUND, str(exc)) from None
