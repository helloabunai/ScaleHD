"""Alembic migrations when the server starts."""

import sqlite3
from pathlib import Path

import pytest
import sqlalchemy as sa
from alembic.autogenerate import compare_metadata
from alembic.operations import Operations
from alembic.runtime.migration import MigrationContext
from fastapi.testclient import TestClient
from scalehd_server.app import create_app
from scalehd_server.config import ServerSettings
from scalehd_server.db import OutdatedDatabaseError, make_engine
from scalehd_server.migrate import head, migrate, migration_connection
from scalehd_server.models import Base

#before migrating
FROZEN = Path(__file__).parent / "data" / "before-migrations.sql"


def _url(path: Path) -> str:
    return f"sqlite:///{path}"


def _load(path: Path) -> None:
    connection = sqlite3.connect(path)
    connection.executescript(FROZEN.read_text())
    connection.close()


def _query(path: Path, sql: str) -> list[tuple[object, ...]]:
    connection = sqlite3.connect(path)
    try:
        return connection.execute(sql).fetchall()
    finally:
        connection.close()


def _tables(path: Path) -> set[str]:
    return {name for (name,) in _query(path, "SELECT name FROM sqlite_master WHERE type='table'")}


def _revision(path: Path) -> str | None:
    if "alembic_version" not in _tables(path):
        return None
    return str(_query(path, "SELECT version_num FROM alembic_version")[0][0])


def _facts(path: Path) -> dict[str, list[tuple[object, ...]]]:
    """What the frozen database contains."""
    return {
        "users": _query(path, "SELECT username, is_admin FROM users ORDER BY id"),
        "login_sessions": _query(path, "SELECT user_id, token_hash FROM login_sessions"),
        "tags": _query(path, "SELECT name, created_by_id FROM tags ORDER BY id"),
        "jobs": _query(path, "SELECT owner_id, name, status FROM jobs ORDER BY id"),
        "job_tags": _query(path, "SELECT job_id, tag_id FROM job_tags ORDER BY job_id, tag_id"),
        "samples": _query(
            path, "SELECT job_id, name, status, genotype, error FROM samples ORDER BY id"
        ),
    }


def test_a_new_database_is_built_by_the_migrations_to_match_the_models(tmp_path: Path) -> None:
    database = tmp_path / "scalehd.db"
    migrate(_url(database))
    assert _revision(database) == head()
    engine = make_engine(_url(database))
    with engine.connect() as connection:
        differences = compare_metadata(MigrationContext.configure(connection), Base.metadata)
    engine.dispose()
    assert differences == []
    assert list(tmp_path.glob("*.before-*")) == []


def test_a_database_from_before_migrations_is_taken_over_with_its_rows(tmp_path: Path) -> None:
    database = tmp_path / "scalehd.db"
    _load(database)
    before = _facts(database)
    migrate(_url(database))
    assert _revision(database) == head()
    assert _facts(database) == before
    assert _query(database, "PRAGMA foreign_key_check") == []
    copy = tmp_path / f"scalehd.db.before-{head()}"
    assert _revision(copy) is None
    assert _facts(copy) == before


def test_a_database_from_before_tags_gets_its_tags_tables(tmp_path: Path) -> None:
    database = tmp_path / "scalehd.db"
    _load(database)
    connection = sqlite3.connect(database)
    connection.executescript("DROP TABLE job_tags; DROP TABLE tags;")
    connection.close()
    migrate(_url(database))
    assert {"tags", "job_tags"} <= _tables(database)
    assert _query(database, "SELECT name FROM jobs ORDER BY id") == [
        ("run-01",),
        ("cancelled run",),
    ]


def test_a_database_missing_columns_is_refused_and_left_as_it_was(tmp_path: Path) -> None:
    database = tmp_path / "scalehd.db"
    connection = sqlite3.connect(database)
    connection.execute("CREATE TABLE users (id INTEGER PRIMARY KEY, username TEXT)")
    connection.commit()
    connection.close()
    with pytest.raises(OutdatedDatabaseError, match=r"scalehd\.db is from an older ScaleHD"):
        migrate(_url(database))
    assert _tables(database) == {"users"}


def test_a_database_already_up_to_date_is_not_copied_again(tmp_path: Path) -> None:
    database = tmp_path / "scalehd.db"
    _load(database)
    migrate(_url(database))
    copies = sorted(tmp_path.glob("*.before-*"))
    migrate(_url(database))
    assert sorted(tmp_path.glob("*.before-*")) == copies


def test_rebuilding_a_table_keeps_the_rows_that_point_at_it(tmp_path: Path) -> None:
    database = tmp_path / "scalehd.db"
    _load(database)
    migrate(_url(database))
    with migration_connection(_url(database)) as connection, connection.begin():
        operations = Operations(MigrationContext.configure(connection))
        with operations.batch_alter_table("jobs", recreate="always") as batch:
            batch.add_column(sa.Column("note", sa.String(), nullable=True))
    assert len(_query(database, "SELECT * FROM job_tags")) == 2
    assert len(_query(database, "SELECT * FROM samples")) == 3
    assert _query(database, "PRAGMA foreign_key_check") == []


def test_the_server_brings_its_database_up_to_date_on_start(tmp_path: Path) -> None:
    settings = ServerSettings(database_dir=tmp_path / "data", workspace=tmp_path / "ws", workers=1)
    with TestClient(create_app(settings)) as client:
        assert client.get("/api/health").status_code == 200
    assert _revision(tmp_path / "data" / "scalehd.db") == head()
