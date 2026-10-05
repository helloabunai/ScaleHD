from collections.abc import Callable, Iterator
from concurrent.futures import Executor, Future
from pathlib import Path
from typing import Any

import pytest
from fastapi.testclient import TestClient
from scalehd_server.app import create_app
from scalehd_server.config import ServerSettings
from scalehd_server.db import make_engine
from scalehd_server.models import Base
from sqlalchemy.orm import Session, sessionmaker


@pytest.fixture
def client(tmp_path: Path) -> Iterator[TestClient]:
    with TestClient(
        create_app(
            ServerSettings(
                database_dir=tmp_path / "data", workspace=tmp_path / "workspace", workers=1
            )
        )
    ) as client:
        yield client


@pytest.fixture
def logged_in(client: TestClient) -> TestClient:
    """The client, logged in as autotest-user, the first account and so the admin."""
    response = client.post(
        "/api/auth/register", json={"username": "autotest-user", "password": "correct horse"}
    )
    assert response.status_code == 201
    return client


@pytest.fixture
def sessions(tmp_path: Path) -> Iterator[sessionmaker[Session]]:
    """Sessions on a fresh SQLite database with every table."""
    engine = make_engine(f"sqlite:///{tmp_path / 'test.db'}")
    Base.metadata.create_all(engine)
    yield sessionmaker(engine)
    engine.dispose()


class StandInPool(Executor):
    """Runs nothing. Each submitted sample waits for the test to finish."""

    def __init__(self) -> None:
        self.submitted: list[tuple[Any, Future[Any]]] = []
        self.shut_down = False

    def submit(self, fn: Any, /, *args: Any, **kwargs: Any) -> Future[Any]:
        future: Future[Any] = Future()
        self.submitted.append((args[0], future))
        return future

    def shutdown(self, wait: bool = True, *, cancel_futures: bool = False) -> None:
        self.shut_down = True


class StandInPools:
    """A pool factory for JobRunner that keeps every pool it makes."""

    def __init__(self) -> None:
        self.made: list[StandInPool] = []

    def __call__(self, workers: int) -> StandInPool:
        self.made.append(StandInPool())
        return self.made[-1]

    @property
    def submitted(self) -> list[tuple[Any, Future[Any]]]:
        return [entry for pool in self.made for entry in pool.submitted]


@pytest.fixture
def pools() -> StandInPools:
    return StandInPools()


@pytest.fixture
def quiet_client(tmp_path: Path, pools: StandInPools) -> Iterator[TestClient]:
    """Logged in as autotest-user. jobs go to standby pools, so no real runs."""
    settings = ServerSettings(
        database_dir=tmp_path / "data", workspace=tmp_path / "workspace", workers=1
    )
    with TestClient(create_app(settings, executor_factory=pools)) as client:
        response = client.post(
            "/api/auth/register", json={"username": "autotest-user", "password": "correct horse"}
        )
        assert response.status_code == 201
        yield client


@pytest.fixture
def data_root(tmp_path: Path) -> Path:
    """A faked data folder with one run of samples.
    Sample a = paired, sample b with R1 only (empty files)."""
    root = tmp_path / "data root"
    run = root / "run-01"
    run.mkdir(parents=True)
    for name in ("a_R1.fastq.gz", "a_R2.fastq.gz", "b_R1.fastq.gz"):
        (run / name).write_bytes(b"")
    return root


@pytest.fixture
def queued(tmp_path: Path, data_root: Path, pools: StandInPools) -> Iterator[TestClient]:
    """Logged in as autotest-user (the admin), with ``data_root`` as the data folder.
    Doesn't actually run just mocked.
    """
    settings = ServerSettings(
        database_dir=tmp_path / "db",
        workspace=tmp_path / "workspace",
        data_root=data_root,
        workers=1,
    )
    with TestClient(create_app(settings, executor_factory=pools)) as client:
        body = {"username": "autotest-user", "password": "correct horse"}
        assert client.post("/api/auth/register", json=body).status_code == 201
        yield client


@pytest.fixture
def login_as() -> Callable[[TestClient, str], None]:
    """Switch a client to another user, registering them the first time."""

    def switch(client: TestClient, username: str) -> None:
        client.cookies.clear()
        body = {"username": username, "password": "correct horse"}
        if client.post("/api/auth/login", json=body).status_code != 200:
            assert client.post("/api/auth/register", json=body).status_code == 201

    return switch
