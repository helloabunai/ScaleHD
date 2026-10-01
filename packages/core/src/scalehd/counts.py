"""Per-sample tallies of molecule repeat structures."""

from __future__ import annotations

import json
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from .pairs import DEFAULT_PREFER_R1, DiscordancePolicy, join_read_base_pairings
from .parse import ReadParse, RepeatParser
from .seqio import read_fastq, read_pairs, reverse_complement
from .structure import FIELDS, AlleleStructure, FieldStatus, Observation

SCHEMA = "scalehd.counts/1"


@dataclass
class SampleCounts:
    """Molecule structures for one sample.

    ``complete`` holds molecules whose every repeat tract was seen in full; ``partial`` holds
    the rest (typically long alleles whose reads end inside the CAG tract).
    ``dropped`` counts molecules discarded because their read base pairings disagreed.
    """

    complete: Counter[AlleleStructure] = field(default_factory=Counter)
    partial: Counter[Observation] = field(default_factory=Counter)
    read_outcomes: Counter[str] = field(default_factory=Counter)
    discordant: Counter[str] = field(default_factory=Counter)
    molecules: int = 0
    unusable: int = 0
    dropped: int = 0

    def add(self, observation: Observation | None, *, dropped: bool = False) -> None:
        self.molecules += 1
        if dropped:
            self.dropped += 1
        elif observation is None or observation.is_empty:
            self.unusable += 1
        elif observation.is_complete:
            self.complete[observation.structure()] += 1
        else:
            self.partial[observation] += 1

    def top(self, n: int = 10) -> list[tuple[AlleleStructure, int]]:
        return self.complete.most_common(n)

    def cag_ccg_matrix(self, max_cag: int = 200, max_ccg: int = 20) -> np.ndarray:
        """Complete molecules as a (CCG, CAG) count matrix, 1-based along both axes.

        Same layout as the ScaleHD 1.x forward distribution. Strucutres longer than
        what could be read are dropped (for now?)
        """
        matrix = np.zeros((max_ccg, max_cag), dtype=np.int64)
        for structure, n in self.complete.items():
            if 1 <= structure.cag <= max_cag and 1 <= structure.ccg <= max_ccg:
                matrix[structure.ccg - 1, structure.cag - 1] += n
        return matrix

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "molecules": self.molecules,
            "unusable": self.unusable,
            "dropped": self.dropped,
            "complete": [
                {"structure": s.label, "count": n} for s, n in self.complete.most_common()
            ],
            "partial": [
                {
                    "observation": o.label,
                    "counts": list(o.counts),
                    "status": [s.name.lower() for s in o.status],
                    "count": n,
                }
                for o, n in self.partial.most_common()
            ],
            "read_outcomes": dict(sorted(self.read_outcomes.items())),
            "discordant": dict(sorted(self.discordant.items())),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> SampleCounts:
        if data.get("schema") != SCHEMA:
            raise ValueError(f"unsupported counts schema {data.get('schema')!r}")
        partial: Counter[Observation] = Counter()
        for row in data["partial"]:
            status = tuple(FieldStatus[s.upper()] for s in row["status"])
            partial[Observation(tuple(row["counts"]), status)] = row["count"]
        return cls(
            complete=Counter(
                {AlleleStructure.from_label(r["structure"]): r["count"] for r in data["complete"]}
            ),
            partial=partial,
            read_outcomes=Counter(data["read_outcomes"]),
            discordant=Counter(data["discordant"]),
            molecules=data["molecules"],
            unusable=data["unusable"],
            dropped=data["dropped"],
        )

    def write_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2) + "\n")

    @classmethod
    def read_json(cls, path: str | Path) -> SampleCounts:
        return cls.from_dict(json.loads(Path(path).read_text()))


def _usable(parse: ReadParse) -> Observation | None:
    return parse.observation if parse.usable else None


def count_reads(
    reads: Iterable[tuple[str, str | None]],
    parser: RepeatParser | None = None,
    policy: DiscordancePolicy = DiscordancePolicy.DROP,
    prefer_r1: tuple[bool, ...] = DEFAULT_PREFER_R1,
) -> SampleCounts:
    """Tally ``(r1, r2)`` sequence pairs as sequenced; pass ``r2=None`` for single-end."""
    parser = parser or RepeatParser()
    counts = SampleCounts()
    for r1, r2 in reads:
        first = parser.parse(r1)
        counts.read_outcomes[f"r1.{first.outcome}"] += 1
        if r2 is None:
            counts.add(_usable(first))
            continue
        second = parser.parse(reverse_complement(r2))
        counts.read_outcomes[f"r2.{second.outcome}"] += 1
        joined = join_read_base_pairings(
            _usable(first), _usable(second), policy=policy, prefer_r1=prefer_r1
        )
        for name, flag in zip(FIELDS, joined.discordant, strict=True):
            if flag:
                counts.discordant[name] += 1
        dropped = joined.is_discordant and policy is DiscordancePolicy.DROP
        counts.add(joined.observation, dropped=dropped)
    return counts


def count_fastq(
    r1: str | Path,
    r2: str | Path | None = None,
    parser: RepeatParser | None = None,
    policy: DiscordancePolicy = DiscordancePolicy.DROP,
    prefer_r1: tuple[bool, ...] = DEFAULT_PREFER_R1,
) -> SampleCounts:
    if r2 is None:
        pairs: Iterable[tuple[str, str | None]] = ((rec.sequence, None) for rec in read_fastq(r1))
    else:
        pairs = ((a.sequence, b.sequence) for a, b in read_pairs(r1, r2))
    return count_reads(pairs, parser, policy, prefer_r1)
