#!/usr/bin/env bash
# Replace legacy/ with ScaleHD 1.x as maintained upstream (MoncktonLab), at one commit, unmodified.
#
#   tools/update-legacy.sh <commit>
#
# Writes legacy/UPSTREAM (where it came from) and legacy/UPSTREAM.sha256 (a checksum per file,
# checked by apps/server/tests/test_legacy_snapshot.py so nothing in legacy/ is edited by hand).
set -euo pipefail

UPSTREAM=https://github.com/MoncktonLab/ScaleHD
COMMIT=${1:?usage: tools/update-legacy.sh <commit>}
REPO_ROOT=$(cd "$(dirname "$0")/.." && pwd)
LEGACY="$REPO_ROOT/legacy"

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
git clone --quiet "$UPSTREAM" "$work/upstream"
full=$(git -C "$work/upstream" rev-parse --verify "$COMMIT^{commit}")
date=$(git -C "$work/upstream" log -1 --format=%cs "$full")

# Upstream keeps the package under src/. Here it sits at legacy/ directly, as before. Its
# tests (test/, the lab's ContinuousTesting/) come too, where that commit has them.
paths=(src pytest.ini README.md)
for optional in test ContinuousTesting; do
    if git -C "$work/upstream" cat-file -e "$full:$optional" 2>/dev/null; then
        paths+=("$optional")
    fi
done
mkdir "$work/snapshot"
git -C "$work/upstream" archive "$full" "${paths[@]}" | tar -x -C "$work/snapshot"
rm -rf "$LEGACY"
mkdir "$LEGACY"
cp -a "$work/snapshot/src/." "$LEGACY/"
for path in "${paths[@]:1}"; do
    if [[ "$path" == README.md ]]; then
        cp -a "$work/snapshot/README.md" "$LEGACY/README.repo-root.md"
    else
        cp -a "$work/snapshot/$path" "$LEGACY/"
    fi
done

printf 'url=%s\ncommit=%s\ndate=%s\n' "$UPSTREAM" "$full" "$date" > "$LEGACY/UPSTREAM"
(cd "$LEGACY" && find . -type f ! -name UPSTREAM ! -name UPSTREAM.sha256 -print0 \
    | LC_ALL=C sort -z | xargs -0 sha256sum) > "$LEGACY/UPSTREAM.sha256"

echo "legacy/ is now $UPSTREAM at ${full:0:7} ($date)"
