"""Database engine and per-request sessions (SQLAlchemy 2)."""

from __future__ import annotations

import sqlite3
from collections.abc import Iterator
from datetime import UTC, datetime
from typing import Annotated, Any, ClassVar

from fastapi import Depends, Request
from sqlalchemy import DateTime, Dialect, Engine, TypeDecorator, create_engine, event, inspect
from sqlalchemy.orm import DeclarativeBase, Session


class UTCDateTime(TypeDecorator[datetime]):
    """Datetimes stored as UTC and read back timezone-aware, SQLite included."""

    impl = DateTime
    cache_ok = True

    def process_bind_param(self, value: datetime | None, dialect: Dialect) -> datetime | None:
        return None if value is None else value.astimezone(UTC).replace(tzinfo=None)

    def process_result_value(self, value: Any, dialect: Dialect) -> datetime | None:
        return None if value is None else value.replace(tzinfo=UTC)


class Base(DeclarativeBase):
    type_annotation_map: ClassVar[dict[Any, Any]] = {datetime: UTCDateTime}


def make_engine(url: str) -> Engine:
    if not url.startswith("sqlite"):
        return create_engine(url)
    engine = create_engine(url, connect_args={"check_same_thread": False})

    @event.listens_for(engine, "connect")
    def _pragmas(connection: sqlite3.Connection, _record: object) -> None:
        # SQLite ignores foreign keys unless asked. WAL lets the web interface read
        # while the job runner writes.
        cursor = connection.cursor()
        cursor.execute("PRAGMA foreign_keys=ON")
        cursor.execute("PRAGMA journal_mode=WAL")
        # Request handlers and the job runner's callback thread both write: wait up
        # to 5 s for the other instead of failing.
        cursor.execute("PRAGMA busy_timeout=5000")
        cursor.close()

    return engine


class OutdatedDatabaseError(RuntimeError):
    """The database was made by an older ScaleHD and lacks columns this version needs."""


def check_schema(engine: Engine) -> None:
    """Stop with a clear message if existing tables lack columns the models have.

    create_all only creates missing tables; it can't add columns to old ones.
    TODO: Alembic migrations before releasing.
    """
    inspector = inspect(engine)
    for table in Base.metadata.sorted_tables:
        if not inspector.has_table(table.name):
            continue
        existing = {column["name"] for column in inspector.get_columns(table.name)}
        if any(column.name not in existing for column in table.columns):
            where = (
                engine.url.database
                if engine.url.get_backend_name() == "sqlite"
                else engine.url.render_as_string(hide_password=True)
            )
            raise OutdatedDatabaseError(
                f"database at {where} is from an older ScaleHD version: "
                "delete it and restart (all accounts will be lost)"
            )


def get_session(request: Request) -> Iterator[Session]:
    with request.app.state.sessions() as session:
        yield session


DbSession = Annotated[Session, Depends(get_session)]
