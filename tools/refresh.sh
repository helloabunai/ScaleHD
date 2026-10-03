#!/usr/bin/env bash
# Rebuild the ScaleHD Docker stack from scratch.
#
# Stops the containers, deletes the database volume, then rebuilds and starts
# the stack. Completely fresh server. Pass -y to skip the confirmation prompt.
#
#   tools/refresh.sh        # asks before deleting
#   tools/refresh.sh -y     # no prompt

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VOLUME="scalehd_scalehd-data"

cd "$REPO_ROOT"

cat >&2 <<WARN

  WARNING: this deletes everything stored on the ScaleHD Docker server.

  Removes the database volume "$VOLUME"!! all user accounts, jobs, settings
  and run history are gone and cannot be recovered.

  Only the Docker volume is deleted. Files under SCALEHD_WORKSPACE and
  SCALEHD_DATA_ROOT on the host are left alone.

WARN

if [[ "${1:-}" != "-y" && "${1:-}" != "--yes" ]]; then
    read -r -p "Type 'yes' to continue: " answer
    if [[ "$answer" != "yes" ]]; then
        echo "Aborted." >&2
        exit 1
    fi
fi

echo "==> docker compose down"
docker compose down

echo "==> docker volume rm $VOLUME"
if docker volume inspect "$VOLUME" >/dev/null 2>&1; then
    docker volume rm "$VOLUME"
else
    echo "    (volume not found, nothing to remove)"
fi

echo "==> docker compose up --build -d"
docker compose up --build -d

echo
echo "Done. ScaleHD is starting at http://localhost:8000"
