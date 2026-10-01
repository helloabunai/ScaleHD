# scalehd-server

The HTTP API, job runner and database behind the ScaleHD web interface, built on the
`scalehd` core. A skeleton for now: every route is declared with its request and
response types, but only `/api/health` works. The rest answer 501 "not implemented yet".

```sh
uv run scalehd-server --reload     # http://127.0.0.1:8000, API docs at /api/docs
```

Settings are `SCALEHD_*` environment variables (`config.py`): `DATA_DIR`,
`DATABASE_URL`, `INPUT_DIR`, `WEB_DIR`, `WORKERS`, `SECRET_KEY`, `ALLOW_REGISTRATION`.

| module | what | state |
|---|---|---|
| `app.py` | app factory: database, runner, `/api` routes, frontend at `/` | works |
| `config.py` | settings from the environment | works |
| `db.py`, `models.py` | SQLAlchemy engine, sessions; users, jobs, samples tables | works, tables created at startup (no migrations yet) |
| `schemas.py` | API request and response bodies, job settings → `CallerSettings` | works |
| `runner.py` | `run_sample` (count and call one sample); `JobRunner` process pool | `run_sample` works, queueing is a stub |
| `auth.py` | password hashing, session cookie, current user | stub, plan in its docstring |
| `routes/` | `health`, `accounts`, `inputs`, `jobs`, `settings` | only `health` works |
