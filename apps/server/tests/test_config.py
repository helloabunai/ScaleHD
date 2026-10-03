"""Settings come from SCALEHD_* variables. Host machine folders are made
absolute in docker server."""

import sqlite3
from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from scalehd_server.app import create_app
from scalehd_server.config import ServerSettings
from scalehd_server.db import OutdatedDatabaseError


def test_settings_come_from_scalehd_variables(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setenv("SCALEHD_DATABASE_DIR", str(tmp_path / "db"))
    monkeypatch.setenv("SCALEHD_DATA_ROOT", str(tmp_path / "runs"))
    monkeypatch.setenv("SCALEHD_WORKSPACE", str(tmp_path / "results"))
    settings = ServerSettings()
    assert settings.database_dir == tmp_path / "db"
    assert settings.data_root == (tmp_path / "runs").resolve()
    assert settings.workspace == (tmp_path / "results").resolve()


def test_no_data_root_by_default(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("SCALEHD_DATA_ROOT", raising=False)
    assert ServerSettings().data_root is None


def test_a_relative_workspace_is_made_absolute(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.chdir(tmp_path)
    workspace = ServerSettings(workspace=Path("results")).workspace
    assert workspace.is_absolute()
    assert workspace == (tmp_path / "results").resolve()


def test_database_from_an_older_version_stops_startup(tmp_path: Path) -> None:
    database_dir = tmp_path / "data"
    database_dir.mkdir()
    old = sqlite3.connect(database_dir / "scalehd.db")
    old.execute("CREATE TABLE users (id INTEGER PRIMARY KEY, username TEXT)")
    old.commit()
    old.close()
    settings = ServerSettings(database_dir=database_dir, workspace=tmp_path / "ws", workers=1)
    with (
        pytest.raises(OutdatedDatabaseError, match=r"scalehd\.db is from an older ScaleHD version"),
        TestClient(create_app(settings)),
    ):
        pass
