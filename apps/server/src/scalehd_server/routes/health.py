"""Is the server up, and which versions is it running."""

import platform
import sqlite3
from functools import cache
from importlib.metadata import version

import scalehd
from fastapi import APIRouter

from .. import __version__
from ..schemas import Health

router = APIRouter(tags=["health"])

# numpy and scipy are probably scientifically relevant to genotype calls
# so make it available from API. TODO: include numpy/scipy ver in run reports when implemented
LIBRARIES = ("numpy", "scipy", "fastapi")


@cache
def _health() -> Health:
    # Docker's healthcheck asks every 30 seconds.
    return Health(
        status="ok",
        version=__version__,
        core_version=scalehd.__version__,
        python=platform.python_version(),
        platform=f"{platform.system()} {platform.machine()}",
        sqlite=sqlite3.sqlite_version,
        libraries={name: version(name) for name in LIBRARIES},
    )


@router.get("/health")
def health() -> Health:
    return _health()
