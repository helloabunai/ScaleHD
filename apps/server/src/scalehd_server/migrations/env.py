"""Run alembic for migrations. ``migrate.py`` or command line
(``tools/new-migration.sh``) with ``-x url=db_url_here``."""

from __future__ import annotations

from typing import Any, Literal

from alembic import context
from alembic.autogenerate.api import AutogenContext
from scalehd_server import models  # noqa: F401
from scalehd_server.db import Base, UTCDateTime
from scalehd_server.migrate import migration_connection
from sqlalchemy import Connection


def _render_item(type_: str, obj: Any, autogen_context: AutogenContext) -> str | Literal[False]:
    if type_ == "type" and isinstance(obj, UTCDateTime):
        return "sa.DateTime()"
    return False


def _run(connection: Connection) -> None:
    context.configure(
        connection=connection,
        target_metadata=Base.metadata,
        render_as_batch=True,  # SQLite changes a table by rebuilding it
        render_item=_render_item,
    )
    with context.begin_transaction():
        context.run_migrations()
        broken = connection.exec_driver_sql("PRAGMA foreign_key_check").fetchall()
        if broken:
            raise RuntimeError(f"the migration left rows pointing at nothing: {broken}")


if context.is_offline_mode():
    raise SystemExit("ScaleHD's migrations run against a database, not to SQL")
handed = context.config.attributes.get("connection")
if handed is not None:
    _run(handed)
else:
    url = context.get_x_argument(as_dictionary=True).get("url")
    if url is None:
        raise SystemExit("name the database: alembic -x url=sqlite:///path/to.db ...")
    with migration_connection(url) as connection:
        _run(connection)
