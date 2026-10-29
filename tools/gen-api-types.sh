#!/usr/bin/env bash
# Generate API types (apps/web/src/api.gen.ts)
# Run this when changing schema/route in app/server
#
#   tools/gen-api-types.sh [--check]
#
# --check only compares (for wip testing)
#
# openapi-typescript lives in tools/api-types, not apps/web
# it needs the TypeScript 5 API, which apps/web's TS7 doesn't have yet

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOOL="$REPO_ROOT/tools/api-types"
OUT="apps/web/src/api.gen.ts"
check=false
if [[ "${1:-}" == "--check" ]]; then
    check=true
fi

cd "$REPO_ROOT"

# Install, or reinstall when package-lock.json changed since the last install (as in
# tools/dev.sh for apps/web).
if [[ ! -f "$TOOL/node_modules/.package-lock.json" ||
    "$TOOL/package-lock.json" -nt "$TOOL/node_modules/.package-lock.json" ]]; then
    echo "==> npm ci (tools/api-types)"
    (cd "$TOOL" && npm ci --no-audit --no-fund)
fi

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

uv run python -m scalehd_server.openapi >"$tmp/openapi.json"
"$TOOL/node_modules/.bin/openapi-typescript" "$tmp/openapi.json" --output "$tmp/api.gen.ts" >/dev/null

if [[ "$check" == true ]]; then
    if ! cmp -s "$tmp/api.gen.ts" "$OUT"; then
        echo "$OUT is out of date: run tools/gen-api-types.sh" >&2
        exit 1
    fi
    echo "$OUT is up to date."
else
    mv "$tmp/api.gen.ts" "$OUT"
    echo "Wrote $OUT"
fi
