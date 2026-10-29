"""The schema at the time of implementing migrations (main @ commit 848b161).

aka first actual DB state going forward

Revision ID: 0001
Revises:
"""

import sqlalchemy as sa
from alembic import op
from scalehd_server.db import OutdatedDatabaseError

revision = "0001"
down_revision = None
branch_labels = None
depends_on = None

_COLUMNS = {
    "users": (
        "id",
        "username",
        "password_hash",
        "is_admin",
        "created_at",
        "default_settings",
        "theme",
    ),
    "login_sessions": ("id", "user_id", "token_hash", "created_at", "expires_at"),
    "tags": ("id", "name", "created_by_id", "created_at"),
    "jobs": (
        "id",
        "owner_id",
        "name",
        "status",
        "settings",
        "created_at",
        "started_at",
        "finished_at",
        "error",
        "demo",
        "output_dir",
    ),
    "job_tags": ("job_id", "tag_id"),
    "samples": (
        "id",
        "job_id",
        "name",
        "r1",
        "r2",
        "simulation",
        "truth",
        "matches_truth",
        "status",
        "genotype",
        "confidence",
        "flags",
        "call",
        "error",
    ),
}


def upgrade() -> None:
    bind = op.get_bind()
    existing = set(sa.inspect(bind).get_table_names())
    for name, make in _MAKE.items():
        if name not in existing:
            make()
            continue
        found = {column["name"] for column in sa.inspect(bind).get_columns(name)}
        if not set(_COLUMNS[name]) <= found:
            raise OutdatedDatabaseError(
                f"database at {bind.engine.url.database} is from an older ScaleHD version: "
                "delete it and restart (all accounts will be lost)"
            )


def downgrade() -> None:
    raise NotImplementedError("migrations only go forwards: restore the copy made before it")


def _users() -> None:
    op.create_table(
        "users",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("username", sa.String(length=64), nullable=False),
        sa.Column("password_hash", sa.String(), nullable=False),
        sa.Column("is_admin", sa.Boolean(), nullable=False),
        sa.Column("created_at", sa.DateTime(), nullable=False),
        sa.Column("default_settings", sa.JSON(), nullable=False),
        sa.Column("theme", sa.Enum("SYSTEM", "LIGHT", "DARK", name="theme"), nullable=False),
        sa.PrimaryKeyConstraint("id", name=op.f("pk_users")),
        sa.UniqueConstraint("username", name=op.f("uq_users_username")),
    )


def _login_sessions() -> None:
    op.create_table(
        "login_sessions",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("user_id", sa.Integer(), nullable=False),
        sa.Column("token_hash", sa.String(length=64), nullable=False),
        sa.Column("created_at", sa.DateTime(), nullable=False),
        sa.Column("expires_at", sa.DateTime(), nullable=False),
        sa.ForeignKeyConstraint(
            ["user_id"],
            ["users.id"],
            name=op.f("fk_login_sessions_user_id_users"),
            ondelete="CASCADE",
        ),
        sa.PrimaryKeyConstraint("id", name=op.f("pk_login_sessions")),
        sa.UniqueConstraint("token_hash", name=op.f("uq_login_sessions_token_hash")),
    )


def _tags() -> None:
    op.create_table(
        "tags",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("name", sa.String(length=15), nullable=False),
        sa.Column("created_by_id", sa.Integer(), nullable=True),
        sa.Column("created_at", sa.DateTime(), nullable=False),
        sa.ForeignKeyConstraint(
            ["created_by_id"], ["users.id"], name=op.f("fk_tags_created_by_id_users")
        ),
        sa.PrimaryKeyConstraint("id", name=op.f("pk_tags")),
        sa.UniqueConstraint("name", name=op.f("uq_tags_name")),
    )


def _jobs() -> None:
    op.create_table(
        "jobs",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("owner_id", sa.Integer(), nullable=False),
        sa.Column("name", sa.String(length=200), nullable=False),
        sa.Column(
            "status",
            sa.Enum(
                "QUEUED",
                "RUNNING",
                "CANCELLING",
                "FINISHED",
                "FAILED",
                "CANCELLED",
                name="jobstatus",
            ),
            nullable=False,
        ),
        sa.Column("settings", sa.JSON(), nullable=False),
        sa.Column("created_at", sa.DateTime(), nullable=False),
        sa.Column("started_at", sa.DateTime(), nullable=True),
        sa.Column("finished_at", sa.DateTime(), nullable=True),
        sa.Column("error", sa.String(), nullable=True),
        sa.Column("demo", sa.Boolean(), nullable=False),
        sa.Column("output_dir", sa.String(), nullable=True),
        sa.ForeignKeyConstraint(["owner_id"], ["users.id"], name=op.f("fk_jobs_owner_id_users")),
        sa.PrimaryKeyConstraint("id", name=op.f("pk_jobs")),
    )


def _job_tags() -> None:
    op.create_table(
        "job_tags",
        sa.Column("job_id", sa.Integer(), nullable=False),
        sa.Column("tag_id", sa.Integer(), nullable=False),
        sa.ForeignKeyConstraint(
            ["job_id"], ["jobs.id"], name=op.f("fk_job_tags_job_id_jobs"), ondelete="CASCADE"
        ),
        sa.ForeignKeyConstraint(
            ["tag_id"], ["tags.id"], name=op.f("fk_job_tags_tag_id_tags"), ondelete="CASCADE"
        ),
        sa.PrimaryKeyConstraint("job_id", "tag_id", name=op.f("pk_job_tags")),
    )


def _samples() -> None:
    op.create_table(
        "samples",
        sa.Column("id", sa.Integer(), nullable=False),
        sa.Column("job_id", sa.Integer(), nullable=False),
        sa.Column("name", sa.String(length=200), nullable=False),
        sa.Column("r1", sa.String(), nullable=True),
        sa.Column("r2", sa.String(), nullable=True),
        sa.Column("simulation", sa.JSON(), nullable=True),
        sa.Column("truth", sa.String(), nullable=True),
        sa.Column("matches_truth", sa.Boolean(), nullable=True),
        sa.Column(
            "status",
            sa.Enum("QUEUED", "RUNNING", "FINISHED", "FAILED", "CANCELLED", name="samplestatus"),
            nullable=False,
        ),
        sa.Column("genotype", sa.String(), nullable=True),
        sa.Column("confidence", sa.Double(), nullable=True),
        sa.Column("flags", sa.JSON(), nullable=False),
        sa.Column("call", sa.JSON(), nullable=True),
        sa.Column("error", sa.String(), nullable=True),
        sa.ForeignKeyConstraint(["job_id"], ["jobs.id"], name=op.f("fk_samples_job_id_jobs")),
        sa.PrimaryKeyConstraint("id", name=op.f("pk_samples")),
    )


_MAKE = {
    "users": _users,
    "login_sessions": _login_sessions,
    "tags": _tags,
    "jobs": _jobs,
    "job_tags": _job_tags,
    "samples": _samples,
}
