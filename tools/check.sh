#!/usr/bin/env bash
# Checks for CI
# --quick arg for git hook
# no arg for full incl slow pytests
#
#   tools/check.sh [--quick]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
quick=false
if [[ "${1:-}" == "--quick" ]]; then
    quick=true
fi

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

if [[ "$quick" == true ]]; then
    echo
    echo "Quick checks passed (tests not run: tools/check.sh runs them)."
    exit 0
fi

echo "==> pytest"
uv run pytest -q

echo
echo "All checks passed."
