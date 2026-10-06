#!/usr/bin/env bash
# Write the next database migration from changes to the server's models e.g.
#
#   tools/new-migration.sh "add a notes column to jobs"
#
# Builds a scratch database from the migrations so far, compares the models with it, and
# writes what changed as the next numbered migration in
# apps/server/src/scalehd_server/migrations/versions/.
# Tests check the migrations still build the right model schema.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VERSIONS="$REPO_ROOT/apps/server/src/scalehd_server/migrations/versions"
message="${1:?usage: tools/new-migration.sh \"what changed\"}"

cd "$REPO_ROOT"
last=$(find "$VERSIONS" -maxdepth 1 -name '[0-9][0-9][0-9][0-9]_*.py' -printf '%f\n' | sort | tail -1)
next=$(printf '%04d' $((10#${last:0:4} + 1)))

scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT
uv run python -c "from scalehd_server.migrate import migrate; migrate('sqlite:///$scratch/latest.db')"
uv run alembic -c apps/server/alembic.ini -x "url=sqlite:///$scratch/latest.db" \
    revision --autogenerate --rev-id "$next" -m "$message"

written=$(find "$VERSIONS" -maxdepth 1 -name "${next}_*.py")
uv run ruff format --quiet "$written"
uv run ruff check --quiet --fix "$written" || true
echo "==> $written: read it, then run the tests"
