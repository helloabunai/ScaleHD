"""The FastAPI application: the HTTP API under ``/api``, the web frontend at ``/``."""

from __future__ import annotations

from collections.abc import AsyncIterator
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from sqlalchemy.orm import sessionmaker
from starlette.exceptions import HTTPException
from starlette.responses import Response
from starlette.types import Scope

from . import __version__, models
from .config import ServerSettings
from .db import OutdatedDatabaseError, check_schema, make_engine
from .routes import api
from .runner import ExecutorFactory, JobRunner, process_pool


class SinglePageApp(StaticFiles):
    """Static files, with ``index.html`` for any other path so frontend routes load."""

    async def get_response(self, path: str, scope: Scope) -> Response:
        try:
            return await super().get_response(path, scope)
        except HTTPException as exc:
            if exc.status_code != 404 or path.split("/")[0] == "api":
                raise
            return await super().get_response("index.html", scope)


def create_app(
    settings: ServerSettings | None = None, *, executor_factory: ExecutorFactory = process_pool
) -> FastAPI:
    settings = settings or ServerSettings()

    @asynccontextmanager
    async def lifespan(app: FastAPI) -> AsyncIterator[None]:
        settings.database_dir.mkdir(parents=True, exist_ok=True)
        engine = make_engine(settings.database)
        try:
            check_schema(engine)
        except OutdatedDatabaseError:
            engine.dispose()
            raise
        # TODO: Alembic migrations once the tables settle.
        models.Base.metadata.create_all(engine)
        sessions = sessionmaker(engine)
        runner = JobRunner(sessions, settings.workers, executor_factory)
        runner.start()
        app.state.settings = settings
        app.state.sessions = sessions
        app.state.runner = runner
        try:
            yield
        finally:
            runner.shutdown()
            engine.dispose()

    app = FastAPI(
        title="ScaleHD",
        version=__version__,
        lifespan=lifespan,
        docs_url="/api/docs",
        redoc_url=None,
        openapi_url="/api/openapi.json",
    )
    app.include_router(api)
    if settings.web_dir is not None:
        app.mount("/", SinglePageApp(directory=settings.web_dir, html=True), name="web")
    return app
