"""The application starts, creates its database, and serves the API and frontend."""

import platform
import sqlite3
from pathlib import Path

import numpy
import scalehd
from fastapi.testclient import TestClient
from scalehd_server import __version__
from scalehd_server.app import create_app
from scalehd_server.config import ServerSettings


def test_health(client: TestClient) -> None:
    response = client.get("/api/health")
    assert response.status_code == 200
    health = response.json()
    assert health["status"] == "ok"
    assert health["version"] == __version__
    assert health["core_version"] == scalehd.__version__


def test_health_says_what_the_server_runs_on(client: TestClient) -> None:
    health = client.get("/api/health").json()
    assert health["python"] == platform.python_version()
    assert health["platform"] == f"{platform.system()} {platform.machine()}"
    assert health["sqlite"] == sqlite3.sqlite_version
    libraries = health["libraries"]
    assert list(libraries) == ["numpy", "scipy", "fastapi"]
    assert libraries["numpy"] == numpy.__version__


def test_database_is_created_in_the_database_directory(client: TestClient, tmp_path: Path) -> None:
    assert (tmp_path / "data" / "scalehd.db").exists()


def test_planned_routes_are_declared(client: TestClient) -> None:
    paths = set(client.get("/api/openapi.json").json()["paths"])
    assert {
        "/api/auth/login",
        "/api/auth/me",
        "/api/inputs",
        "/api/jobs",
        "/api/jobs/demo",
        "/api/jobs/{job_id}",
        "/api/jobs/{job_id}/report",
        "/api/settings",
    } <= paths


def test_routes_need_a_login(client: TestClient) -> None:
    assert client.get("/api/jobs").status_code == 401


def test_routes_not_built_yet_say_so(logged_in: TestClient) -> None:
    response = logged_in.get("/api/jobs/1/report")
    assert response.status_code == 501
    assert "not implemented" in response.json()["detail"]


def test_serves_the_frontend_with_client_side_routes(tmp_path: Path) -> None:
    web = tmp_path / "web"
    web.mkdir()
    (web / "index.html").write_text("<p>ScaleHD</p>")
    settings = ServerSettings(
        database_dir=tmp_path / "data", workspace=tmp_path / "ws", web_dir=web, workers=1
    )
    with TestClient(create_app(settings)) as client:
        assert client.get("/").text == "<p>ScaleHD</p>"
        assert client.get("/jobs/3").text == "<p>ScaleHD</p>"
        assert client.get("/api/no-such-route").status_code == 404
