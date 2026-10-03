#!/usr/bin/env bash
# Run ScaleHD for development, with hot reload etc
#
# Starts the API and Vite dev server (updates
# the page in place when .tsx, .ts or .css files change).
# DEV port is http://localhost:5173, not the docker compose "main" server at port 8000.
#
# Settings come from .env, as usual but database is data/scalehd.db which will be just a dev
# db i.e not the same as in the container/"production" server.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
API="http://127.0.0.1:8000"

cd "$REPO_ROOT"

# Vite forwards /api to port 8000 (apps/web/vite.config.ts), so the API must get it.
if (exec 3<>/dev/tcp/127.0.0.1/8000) 2>/dev/null; then
    echo "Port 8000 is already in use. If it's the Docker stack: docker compose stop" >&2
    exit 1
fi

if [[ -f .env ]]; then
    echo "==> settings from .env"
    set -a
    # shellcheck disable=SC1091
    source .env
    set +a
fi
# A demo run otherwise takes every core.
export SCALEHD_WORKERS="${SCALEHD_WORKERS:-4}"

# Install, or reinstall when package-lock.json changed since the last install (e.g. a
# pulled commit added a dependency). npm writes node_modules/.package-lock.json on install.
if [[ ! -f apps/web/node_modules/.package-lock.json ||
    apps/web/package-lock.json -nt apps/web/node_modules/.package-lock.json ]]; then
    echo "==> npm ci"
    (cd apps/web && npm ci)
fi

server=""
web=""
cleanup() {
    [[ -n "$web" ]] && kill "$web" 2>/dev/null || true
    [[ -n "$server" ]] && kill "$server" 2>/dev/null || true
    wait 2>/dev/null || true
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

echo "==> API at $API (restarts on Python changes)"
uv run scalehd-server --reload &
server=$!

for _ in $(seq 120); do
    curl -sf "$API/api/health" >/dev/null && break
    if ! kill -0 "$server" 2>/dev/null; then
        echo "The API stopped while starting (see above)." >&2
        exit 1
    fi
    sleep 0.5
done
if ! curl -sf "$API/api/health" >/dev/null; then
    echo "The API didn't answer within a minute." >&2
    exit 1
fi

echo "==> web interface (updates on .tsx, .ts and .css changes)"
(cd apps/web && exec ./node_modules/.bin/vite) &
web=$!

# Stop everything if either one exits
while kill -0 "$server" 2>/dev/null && kill -0 "$web" 2>/dev/null; do
    sleep 1
done
