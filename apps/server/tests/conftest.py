from collections.abc import Iterator
from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from scalehd_server.app import create_app
from scalehd_server.config import ServerSettings


@pytest.fixture
def client(tmp_path: Path) -> Iterator[TestClient]:
    with TestClient(create_app(ServerSettings(data_dir=tmp_path, workers=1))) as client:
        yield client


@pytest.fixture
def logged_in(client: TestClient) -> TestClient:
    """The client, logged in as alice, the first account and so the admin."""
    response = client.post(
        "/api/auth/register", json={"username": "alice", "password": "correct horse"}
    )
    assert response.status_code == 201
    return client
