# scalehd-server

The HTTP API, job runner and database behind the ScaleHD web interface, built on the
`scalehd` core. Mostly a skeleton: every route is declared with its request and
response types, but only `/api/health` and accounts work. The rest answer 501 "not
implemented yet" (or 401 when nobody is logged in).

```sh
uv run scalehd-server --reload     # http://127.0.0.1:8000, API docs at /api/docs
```

Settings are `SCALEHD_*` environment variables (`config.py`): `DATA_DIR`,
`DATABASE_URL`, `INPUT_DIR`, `WEB_DIR`, `WORKERS`, `ALLOW_REGISTRATION`, `SESSION_DAYS`.

### Accounts

- The first account created on a new server is the admin. After that, anyone who can
  reach the server can register unless `SCALEHD_ALLOW_REGISTRATION=false`.
- Usernames are case-insensitive (stored in lower case). Passwords need 8 characters
  or more and are stored as argon2 hashes.
- Logging in sets an HttpOnly, SameSite=Lax cookie holding a random token. The
  database keeps only the token's SHA-256, so logging out and expiry
  (`SCALEHD_SESSION_DAYS`, default 30) take effect on the server.
- Changing password logs out every other browser using that account.
- Routes: `GET /api/auth/registration`, `POST /api/auth/register`, `POST /api/auth/login`,
  `POST /api/auth/logout`, `GET /api/auth/me`, `PUT /api/auth/password`. Other routes
  take the logged-in user from the `CurrentUser` dependency in `auth.py`.

| module | what | state |
|---|---|---|
| `app.py` | app factory: database, runner, `/api` routes, frontend at `/` | works |
| `config.py` | settings from the environment | works |
| `db.py`, `models.py` | SQLAlchemy engine, sessions; users, login sessions, jobs, samples tables | works, tables created at startup (no migrations yet) |
| `schemas.py` | API request and response bodies, job settings → `CallerSettings` | works |
| `runner.py` | `run_sample` (count and call one sample); `JobRunner` process pool | `run_sample` works, queueing is a stub |
| `auth.py` | password hashing, session cookie, current user | works |
| `routes/` | `health`, `accounts`, `inputs`, `jobs`, `settings` | `health` and `accounts` work |
