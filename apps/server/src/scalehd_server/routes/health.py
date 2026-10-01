"""Is the server up, and which versions is it running."""

import scalehd
from fastapi import APIRouter

from .. import __version__
from ..schemas import Health

router = APIRouter(tags=["health"])


@router.get("/health")
def health() -> Health:
    return Health(status="ok", version=__version__, core_version=scalehd.__version__)
