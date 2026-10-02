"""Server settings, read from ``SCALEHD_*`` environment variables."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Annotated

from fastapi import Depends, Request
from pydantic import Field, field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict


class ServerSettings(BaseSettings):
    model_config = SettingsConfigDict(env_prefix="SCALEHD_")

    # The SQLite database, unless database_url says otherwise. In Docker it lives in
    # its own volume, apart from the workspace, so nobody browsing results over a share
    # can delete it.
    database_dir: Path = Path("data")
    # SQLAlchemy URL. Unset means SQLite in database_dir, which is plenty for one machine.
    database_url: str | None = None
    # Sequencing data users can browse, read-only. In Docker it is mounted at the same
    # path as on the host, so paths in the web interface match the host.
    data_root: Path | None = None
    # Results: a folder per user, a subfolder per job. Mounted read-write, at the same
    # path as on the host.
    workspace: Path = Path("workspace")
    # The built web frontend (apps/web/dist) to serve at /. Unset in development, where
    # Vite serves the frontend and forwards /api here.
    web_dir: Path | None = None
    # Samples processed at once, one process each.
    workers: int = Field(default_factory=lambda: os.process_cpu_count() or 1, ge=1)
    # Whether anyone who can reach the server may create an account. The first
    # account (the admin) can always be created.
    allow_registration: bool = True
    # How long a login lasts.
    session_days: float = Field(default=30, gt=0)

    @field_validator("data_root", "workspace")
    @classmethod
    def _absolute(cls, folder: Path | None) -> Path | None:
        return None if folder is None else folder.expanduser().resolve()

    @property
    def database(self) -> str:
        return self.database_url or f"sqlite:///{self.database_dir / 'scalehd.db'}"


def get_settings(request: Request) -> ServerSettings:
    settings: ServerSettings = request.app.state.settings
    return settings


ServerConfig = Annotated[ServerSettings, Depends(get_settings)]
