#!/usr/bin/env bash
# Pin ScaleHD 1.x, as maintained upstream (MoncktonLab), at another commit.
#
#   tools/update-legacy.sh <commit>
#
# Writes apps/legacy-runner/UPSTREAM (where 1.x comes from: url, commit, date). The Docker
# image fetches that commit unmodified at build time. Then check the new commit with
# tools/legacy-parity.sh
set -euo pipefail

UPSTREAM=https://github.com/MoncktonLab/ScaleHD
COMMIT=${1:?usage: tools/update-legacy.sh <commit>}
REPO_ROOT=$(cd "$(dirname "$0")/.." && pwd)

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
git clone --quiet --filter=blob:none --no-checkout "$UPSTREAM" "$work/upstream"
full=$(git -C "$work/upstream" rev-parse --verify "$COMMIT^{commit}")
date=$(git -C "$work/upstream" log -1 --format=%cs "$full")

printf 'url=%s\ncommit=%s\ndate=%s\n' "$UPSTREAM" "$full" "$date" \
    > "$REPO_ROOT/apps/legacy-runner/UPSTREAM"
echo "ScaleHD 1.x is now pinned to $UPSTREAM at ${full:0:7} ($date)."
echo "Check it with tools/legacy-parity.sh."
