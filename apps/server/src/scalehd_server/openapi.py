"""API info (generated)

    uv run python -m scalehd_server.openapi

tools/gen-api-types.sh runs this, so the server doesn't need to be running.
"""

from __future__ import annotations

import json

from .app import create_app
from .config import ServerSettings


def main() -> None:
    # No frontend: the description only covers the API.
    app = create_app(ServerSettings(web_dir=None))
    print(json.dumps(app.openapi(), indent=2))


if __name__ == "__main__":
    main()
