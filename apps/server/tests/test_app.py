"""The application starts, creates its database, and serves the API and frontend."""

from collections.abc import Iterator
from pathlib import Path

import pytest
import scalehd
from fastapi.testclient import TestClient
from scalehd_server import __version__
from scalehd_server.app import create_app
from scalehd_server.config import ServerSettings


@pytest.fixture
def client(tmp_path: Path) -> Iterator[TestClient]:
    with TestClient(create_app(ServerSettings(data_dir=tmp_path, workers=1))) as client:
        yield client


def test_health(client: TestClient) -> None:
    response = client.get("/api/health")
    assert response.status_code == 200
    assert response.json() == {
        "status": "ok",
        "version": __version__,
        "core_version": scalehd.__version__,
    }


def test_database_is_created_in_the_data_directory(client: TestClient, tmp_path: Path) -> None:
    assert (tmp_path / "scalehd.db").exists()


def test_planned_routes_are_declared(client: TestClient) -> None:
    paths = set(client.get("/api/openapi.json").json()["paths"])
    assert {
        "/api/auth/login",
        "/api/auth/me",
        "/api/inputs",
        "/api/jobs",
        "/api/jobs/{job_id}",
        "/api/jobs/{job_id}/report",
        "/api/settings",
    } <= paths


def test_routes_not_built_yet_say_so(client: TestClient) -> None:
    response = client.get("/api/jobs")
    assert response.status_code == 501
    assert "not implemented" in response.json()["detail"]


def test_serves_the_frontend_with_client_side_routes(tmp_path: Path) -> None:
    web = tmp_path / "web"
    web.mkdir()
    (web / "index.html").write_text("<p>ScaleHD</p>")
    settings = ServerSettings(data_dir=tmp_path / "data", web_dir=web, workers=1)
    with TestClient(create_app(settings)) as client:
        assert client.get("/").text == "<p>ScaleHD</p>"
        assert client.get("/jobs/3").text == "<p>ScaleHD</p>"
        assert client.get("/api/no-such-route").status_code == 404
