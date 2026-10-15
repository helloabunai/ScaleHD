#!/usr/bin/env bash
# Replace legacy/ with ScaleHD 1.x as maintained upstream (MVLS-RSE), at one commit, unmodified.
#
#   tools/update-legacy.sh <commit>
#
# Writes legacy/UPSTREAM (where it came from) and legacy/UPSTREAM.sha256 (a checksum per file,
# checked by apps/server/tests/test_legacy_snapshot.py so nothing in legacy/ is edited by hand).
set -euo pipefail

UPSTREAM=https://github.com/MVLS-RSE/ScaleHD
COMMIT=${1:?usage: tools/update-legacy.sh <commit>}
REPO_ROOT=$(cd "$(dirname "$0")/.." && pwd)
LEGACY="$REPO_ROOT/legacy"

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
git clone --quiet "$UPSTREAM" "$work/upstream"
full=$(git -C "$work/upstream" rev-parse --verify "$COMMIT^{commit}")
date=$(git -C "$work/upstream" log -1 --format=%cs "$full")

# Upstream keeps the package under src/. Here it sits at legacy/ directly, as before.
mkdir "$work/snapshot"
git -C "$work/upstream" archive "$full" src test pytest.ini README.md | tar -x -C "$work/snapshot"
rm -rf "$LEGACY"
mkdir "$LEGACY"
cp -a "$work/snapshot/src/." "$LEGACY/"
cp -a "$work/snapshot/test" "$work/snapshot/pytest.ini" "$LEGACY/"
cp -a "$work/snapshot/README.md" "$LEGACY/README.repo-root.md"

printf 'url=%s\ncommit=%s\ndate=%s\n' "$UPSTREAM" "$full" "$date" > "$LEGACY/UPSTREAM"
(cd "$LEGACY" && find . -type f ! -name UPSTREAM ! -name UPSTREAM.sha256 -print0 \
    | LC_ALL=C sort -z | xargs -0 sha256sum) > "$LEGACY/UPSTREAM.sha256"

echo "legacy/ is now $UPSTREAM at ${full:0:7} ($date)"
