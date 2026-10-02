#!/usr/bin/env bash
# Run every check CI runs: lint, formatting, types, the web build, then the tests.
#
# Stops at the first failure. The git pre-push hook runs this, and it can be run by
# hand before committing. Checks the working tree as it is, uncommitted changes
# included.
#
#   tools/check.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

cd "$REPO_ROOT"

if [[ ! -d apps/web/node_modules ]]; then
    echo "apps/web/node_modules is missing: run (cd apps/web && npm ci) first." >&2
    exit 1
fi

echo "==> ruff check"
uv run ruff check packages apps

echo "==> ruff format --check"
uv run ruff format --check packages apps

echo "==> mypy"
uv run mypy

echo "==> web type-check and build"
(cd apps/web && npm run build)

echo "==> pytest"
uv run pytest -q

echo
echo "All checks passed."
