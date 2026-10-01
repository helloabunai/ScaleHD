"""``scalehd-server``: run the API (and the built frontend, if configured).

Everything else is configured with ``SCALEHD_*`` environment variables, see ``config.py``.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence

import uvicorn

from . import __version__


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="scalehd-server", description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--reload", action="store_true", help="restart on code changes")
    args = parser.parse_args(argv)
    uvicorn.run(
        "scalehd_server.app:create_app",
        factory=True,
        host=args.host,
        port=args.port,
        reload=args.reload,
    )
    return 0
