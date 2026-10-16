#!/usr/bin/env bash
# Checks the legacy runner gets what ScaleHD 1.x gets on its own, in the image.
# Builds the image first unless one is named.
#
#   tools/legacy-parity.sh [image]
set -euo pipefail

REPO_ROOT=$(cd "$(dirname "$0")/.." && pwd)
IMAGE=${1:-}
if [[ -z "$IMAGE" ]]; then
    IMAGE=scalehd:dev
    docker build -t "$IMAGE" "$REPO_ROOT"
fi
docker run --rm -v "$REPO_ROOT/apps/legacy-runner:/parity:ro" "$IMAGE" python /parity/parity.py
