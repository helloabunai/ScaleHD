"""Database engine and per-request sessions (SQLAlchemy 2)."""

from __future__ import annotations

import sqlite3
from collections.abc import Iterator
from datetime import UTC, datetime
from typing import Annotated, Any, ClassVar

from fastapi import Depends, Request
from sqlalchemy import DateTime, Dialect, Engine, TypeDecorator, create_engine, event
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
        cursor.close()

    return engine


def get_session(request: Request) -> Iterator[Session]:
    with request.app.state.sessions() as session:
        yield session


DbSession = Annotated[Session, Depends(get_session)]
