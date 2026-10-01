"""Database engine and per-request sessions (SQLAlchemy 2)."""

from __future__ import annotations

import sqlite3
from collections.abc import Iterator
from typing import Annotated

from fastapi import Depends, Request
from sqlalchemy import Engine, create_engine, event
from sqlalchemy.orm import DeclarativeBase, Session


class Base(DeclarativeBase):
    pass


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
