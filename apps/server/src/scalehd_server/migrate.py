"""alembic migrations at server boot time. only forward direction.
before upgrades make a copy just incase.
to undo failed migrations, stop the server, put the copy back as the db
"""

from __future__ import annotations

import sqlite3
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

from alembic import command
from alembic.config import Config
from alembic.runtime.migration import MigrationContext
from alembic.script import ScriptDirectory
from sqlalchemy import Connection, create_engine, event, inspect
from sqlalchemy.engine import make_url

SCRIPT_LOCATION = "scalehd_server:migrations"


def _config(connection: Connection | None = None) -> Config:
    config = Config()
    config.set_main_option("script_location", SCRIPT_LOCATION)
    if connection is not None:
        config.attributes["connection"] = connection
    return config


def head() -> str:
    """The latest migration's revision."""
    revision = ScriptDirectory.from_config(_config()).get_current_head()
    if revision is None:
        raise RuntimeError("no migrations found")
    return revision


@contextmanager
def migration_connection(url: str) -> Iterator[Connection]:
    """A connection to migrate the database at ``url`` with."""
    engine = create_engine(url, connect_args={"isolation_level": None})

    @event.listens_for(engine, "connect")
    def _pragmas(dbapi_connection: sqlite3.Connection, _record: object) -> None:
        cursor = dbapi_connection.cursor()
        cursor.execute("PRAGMA foreign_keys=OFF")
        cursor.execute("PRAGMA busy_timeout=5000")
        cursor.close()

    @event.listens_for(engine, "begin")
    def _begin(connection: Connection) -> None:
        connection.exec_driver_sql("BEGIN")

    try:
        with engine.connect() as connection:
            yield connection
    finally:
        engine.dispose()


def migrate(url: str) -> None:
    """Bring the database at ``url`` up to date, create if server is new."""
    if make_url(url).get_backend_name() != "sqlite":
        raise ValueError(f"only SQLite databases are supported for now, not {url!r}")
    latest = head()
    with migration_connection(url) as connection:
        current = MigrationContext.configure(connection).get_current_revision()
        behind = current != latest and bool(inspect(connection).get_table_names())
        connection.rollback()
        if behind:
            _copy(url, f"before-{latest}")
        command.upgrade(_config(connection), "head")


def _copy(url: str, suffix: str) -> None:
    """Copy for backups just in case"""
    database = make_url(url).database
    if not database or database == ":memory:":
        return
    path = Path(database)
    copy = path.with_name(f"{path.name}.{suffix}")
    if copy.exists():
        return
    source, target = sqlite3.connect(path), sqlite3.connect(copy)
    try:
        source.backup(target)
    finally:
        target.close()
        source.close()
