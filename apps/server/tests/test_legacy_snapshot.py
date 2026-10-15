"""legacy/ is ScaleHD 1.x exactly as upstream has it.

Legacy genotyping has to give the same results as ScaleHD 1.x run on its own, so its code
is never edited here. tools/update-legacy.sh replaces it from upstream and records a
checksum per file.
"""

import hashlib
from pathlib import Path

LEGACY = Path(__file__).parents[3] / "legacy"
RECORDS = {"UPSTREAM", "UPSTREAM.sha256"}


def _recorded() -> dict[str, str]:
    lines = (LEGACY / "UPSTREAM.sha256").read_text().splitlines()
    return {
        path.removeprefix("./"): digest for digest, path in (line.split("  ", 1) for line in lines)
    }


def test_every_legacy_file_is_upstreams() -> None:
    changed = [
        path
        for path, digest in _recorded().items()
        if hashlib.sha256((LEGACY / path).read_bytes()).hexdigest() != digest
    ]
    assert changed == [], "edited in legacy/, which must stay as upstream has it"


def test_nothing_is_added_to_or_missing_from_legacy() -> None:
    present = {
        str(path.relative_to(LEGACY))
        for path in LEGACY.rglob("*")
        if path.is_file() and "__pycache__" not in path.parts and path.name not in RECORDS
    }
    assert present == set(_recorded())


def test_legacy_says_where_it_came_from() -> None:
    upstream = dict(line.split("=", 1) for line in (LEGACY / "UPSTREAM").read_text().splitlines())
    assert upstream["url"] == "https://github.com/MVLS-RSE/ScaleHD"
    assert len(upstream["commit"]) == 40
