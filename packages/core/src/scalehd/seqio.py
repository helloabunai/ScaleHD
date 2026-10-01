"""FASTQ reading and writing."""

from __future__ import annotations

import gzip
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import IO

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True, slots=True)
class FastqRecord:
    name: str
    sequence: str
    quality: str

    @property
    def pair_id(self) -> str:
        """Read name with any comment and /1 or /2 suffix removed."""
        name = self.name.split(maxsplit=1)[0]
        if name.endswith(("/1", "/2")):
            name = name[:-2]
        return name


class FastqFormatError(ValueError):
    pass


def open_text(path: str | Path, mode: str = "rt") -> IO[str]:
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="ascii", compresslevel=6)  # type: ignore[return-value]
    return open(path, mode, encoding="ascii")


def read_fastq(path: str | Path) -> Iterator[FastqRecord]:
    with open_text(path) as handle:
        while header := handle.readline():
            sequence = handle.readline().rstrip("\n")
            separator = handle.readline()
            quality = handle.readline().rstrip("\n")
            if not header.startswith("@") or not separator.startswith("+"):
                raise FastqFormatError(f"{path}: malformed record at {header.strip()!r}")
            if len(sequence) != len(quality):
                raise FastqFormatError(f"{path}: sequence/quality length differ in {header!r}")
            yield FastqRecord(header[1:].rstrip("\n"), sequence.upper(), quality)


def read_pairs(r1: str | Path, r2: str | Path) -> Iterator[tuple[FastqRecord, FastqRecord]]:
    """Yield read base pairings, failing if the two files fall out of step."""
    first, second = read_fastq(r1), read_fastq(r2)
    for a in first:
        b = next(second, None)
        if b is None:
            raise FastqFormatError(f"{r2} has fewer reads than {r1}")
        if a.pair_id != b.pair_id:
            raise FastqFormatError(
                f"read base pairings out of sync: {a.pair_id!r} vs {b.pair_id!r}"
            )
        yield a, b
    if next(second, None) is not None:
        raise FastqFormatError(f"{r2} has more reads than {r1}")


def write_fastq(handle: IO[str], records: Iterable[FastqRecord]) -> None:
    for record in records:
        handle.write(f"@{record.name}\n{record.sequence}\n+\n{record.quality}\n")
