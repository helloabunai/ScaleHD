"""Server settings, read from ``SCALEHD_*`` environment variables."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Annotated

from fastapi import Depends, Request
from pydantic import Field, SecretStr
from pydantic_settings import BaseSettings, SettingsConfigDict


class ServerSettings(BaseSettings):
    model_config = SettingsConfigDict(env_prefix="SCALEHD_")

    # Job outputs (counts and calls per sample) and, by default, the SQLite database.
    data_dir: Path = Path("data")
    # SQLAlchemy URL. Unset means SQLite in data_dir, which is plenty for one machine.
    database_url: str | None = None
    # FASTQ files users can pick from in the web interface. Mounted read-only in Docker,
    # so sequencing runs never pass through the browser.
    input_dir: Path | None = None
    # The built web frontend (apps/web/dist) to serve at /. Unset in development, where
    # Vite serves the frontend and forwards /api here.
    web_dir: Path | None = None
    # Samples processed at once, one process each.
    workers: int = Field(default_factory=lambda: os.process_cpu_count() or 1, ge=1)
    # Signs session cookies. Needed once accounts exist.
    secret_key: SecretStr | None = None
    # Whether anyone who can reach the server may create an account.
    allow_registration: bool = True

    @property
    def database(self) -> str:
        return self.database_url or f"sqlite:///{self.data_dir / 'scalehd.db'}"


def get_settings(request: Request) -> ServerSettings:
    settings: ServerSettings = request.app.state.settings
    return settings


ServerConfig = Annotated[ServerSettings, Depends(get_settings)]
