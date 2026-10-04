# scalehd-server

The HTTP API, job runner and database behind the ScaleHD web interface, built on the
`scalehd` core. Mostly a skeleton: every route is declared with its request and
response types,  but most things implemented are just placeholders. Anything else without a placeholder answers 501 "not implemented yet" (or 401 when nobody is logged in).

```sh
uv run scalehd-server --reload     # http://127.0.0.1:8000, API docs at /api/docs
```

Settings are `SCALEHD_*` environment variables (`config.py`): `DATABASE_DIR`,
`DATABASE_URL`, `DATA_ROOT`, `WORKSPACE`, `WEB_DIR`, `WORKERS`, `ALLOW_REGISTRATION`,
`SESSION_DAYS`.

### Folders

- `SCALEHD_DATA_ROOT`: sequencing data users can browse, read-only.
- `SCALEHD_WORKSPACE`: results, `<workspace>/<username>/<job id>-<job name>/`.
- `SCALEHD_DATABASE_DIR`: the SQLite database, kept apart from the workspace.

In Docker the first two are mounted at the same path as on the host, so paths in the
web interface and in result files are host paths. A database made by an older version
stops startup with a message naming it: delete it (accounts are lost) and restart.

### Accounts

- The first account created on a new server is the admin. After that, anyone who can
  reach the server can register unless `SCALEHD_ALLOW_REGISTRATION=false`.
- Usernames are case-insensitive (stored in lower case). Passwords need 8 characters
  or more and are stored as argon2 hashes.
- Logging in sets an HttpOnly, SameSite=Lax cookie holding a random token. The
  database keeps only the token's SHA-256, so logging out and expiry
  (`SCALEHD_SESSION_DAYS`, default 30) take effect on the server.
- Changing password logs out every other browser using that account.
- Each account keeps a light/dark choice (`users.theme`: `system`, `light` or `dark`, default `system`)
- Routes: `GET /api/auth/registration`, `POST /api/auth/register`, `POST /api/auth/login`,
  `POST /api/auth/logout`, `GET /api/auth/me`, `PUT /api/auth/password`,
  `PUT /api/auth/theme`. Other routes take the logged-in user from the `CurrentUser`
  dependency in `auth.py`, or an admin from `AdminUser` (403 for anyone else).
- Admins: `GET /api/admin/users` and `PUT /api/admin/users/{id}/admin` for giving
  existing users admin rights.

### Jobs from the data folder

- `GET /api/inputs?folder=` lists directories within `SCALEHD_DATA_ROOT` and samples within
  those directories.
- `POST /api/jobs` queues selected samples as a job for genotyping.
- Also included UI tags for easy assignment of jobs to e.g. projects / data cohorts.
  Tags can be made on the job page by any user by needs an admin to edit/delete them

### Genotyping method

We plan to allow users to pick between the legacy genotyping method used in ScaleHD 1.x and a newer genotyping approach which may or may not be better, or worse! you're welcome!

Each job runs one of two methods (`JobSettings.method`):

- `legacy`: ScaleHD 1.x, aligning reads to a user provided reference library, then use the 1.x genotyper. The default for new users, but not runnable yet: it still has to be
  extracted from `legacy/`, alignment included. Until then, jobs using it are refused when submitted ("legacy genotyping is not available yet").
- `model`: reads the repeat structure straight from each read, then the model-based
  caller. Shown in the web interface as "New (model-based)", marked Beta.

`GET`/`PUT /api/settings` hold each user's default (stored in `users.default_settings`).
A new job uses its own `settings` if given, otherwise the user's default.
`RUNNABLE_METHODS` in `worker.py` lists what jobs can use.

### Simple job demo

Basically written just to test the end to end/API functionality.

Gets some simulated data from the simulator and runs the new genotyping (only one that exists at the moment) model with data reported back to the frontend.

| module | what | state |
|---|---|---|
| `app.py` | app factory: database, runner, `/api` routes, frontend at `/` | works |
| `config.py` | settings from the environment | works |
| `db.py`, `models.py` | SQLAlchemy engine, sessions; users, login sessions, jobs, samples, tags tables | works, tables created at startup (no migrations yet) |
| `schemas.py` | API request and response bodies, job settings → `CallerSettings` | works |
| `runner.py`, `worker.py` | the job runner (queue, worker pool, recording) and one sample's work | works |
| `workspace.py` | job folders and `job.json` in the workspace | works |
| `demo.py` | the demo job's simulated samples | works |
| `inputs.py` | browsing the data folder, pairing FASTQ files, truth files of simulated runs | works |
| `auth.py` | password hashing, session cookie, current user | works |
| `routes/` | `health`, `accounts`, `admin`, `folders`, `inputs`, `jobs`, `settings`, `tags` | all but the cancel/report `jobs` routes work |
