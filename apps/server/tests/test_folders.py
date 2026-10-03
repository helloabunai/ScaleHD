"""Which server folders contain sequencing data and results, for the web interface."""

from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from scalehd_server.app import create_app
from scalehd_server.config import ServerSettings

TEST_USER = {"username": "autotest-user", "password": "correct horse"}


def test_folders_need_a_login(client: TestClient) -> None:
    assert client.get("/api/folders").status_code == 401


def test_folders_show_the_data_root_workspace_and_your_folder(tmp_path: Path) -> None:
    settings = ServerSettings(
        database_dir=tmp_path / "data",
        data_root=tmp_path / "runs",
        workspace=tmp_path / "ws",
        workers=1,
    )
    with TestClient(create_app(settings)) as client:
        client.post("/api/auth/register", json=TEST_USER)
        folders = client.get("/api/folders").json()
    assert folders == {
        "data_root": str((tmp_path / "runs").resolve()),
        "workspace": str((tmp_path / "ws").resolve()),
        "your_folder": str((tmp_path / "ws" / "autotest-user").resolve()),
    }


def test_no_data_root_is_null(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.delenv("SCALEHD_DATA_ROOT", raising=False)
    settings = ServerSettings(database_dir=tmp_path / "data", workspace=tmp_path / "ws", workers=1)
    with TestClient(create_app(settings)) as client:
        client.post("/api/auth/register", json=TEST_USER)
        assert client.get("/api/folders").json()["data_root"] is None


def test_your_folder_is_where_your_jobs_go(quiet_client: TestClient) -> None:
    job = quiet_client.post("/api/jobs/demo").json()
    your_folder = quiet_client.get("/api/folders").json()["your_folder"]
    assert Path(job["output_dir"]).parent == Path(your_folder)
