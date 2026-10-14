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
from .seqio import FastqRecord, Progress, read_fastq, read_pairs, reverse_complement
from .structure import FIELDS, AlleleStructure, FieldStatus, Observation
from .substitutions import ReadEvidence, Site, molecule_substitutions

# Beta schema 1 = initial implementation idea
# Beta schema 2 adds the "unconfirmed" for extremely long repeats (at the end of sequence length)
# Beta schema 3 adds the evidence for SNPs (only atm substitutions and flank coverage)
SCHEMA = "scalehd.counts/3"
READABLE = ("scalehd.counts/1", "scalehd.counts/2", SCHEMA)


@dataclass
class SampleCounts:
    """Molecule structures for one sample.

    ``complete`` holds molecules whose every repeat tract was seen in full; ``partial`` holds
    the rest (typically long alleles whose reads end inside the CAG tract).
    ``dropped`` counts molecules discarded because their read base pairings disagreed.

    The evidence for SNPs. ``substitutions`` counts molecules showing a substitution, by
    observation, site and base, and ``flank_coverage`` counts molecules by observation and
    whether their reads reached the 5' and 3' flanks. Both are empty for counts made
    before SNPs were called (schema 1 and 2).
    """

    complete: Counter[AlleleStructure] = field(default_factory=Counter)
    partial: Counter[Observation] = field(default_factory=Counter)
    read_outcomes: Counter[str] = field(default_factory=Counter)
    discordant: Counter[str] = field(default_factory=Counter)
    molecules: int = 0
    unusable: int = 0
    dropped: int = 0
    substitutions: Counter[tuple[Observation, Site, str]] = field(default_factory=Counter)
    flank_coverage: Counter[tuple[Observation, bool, bool]] = field(default_factory=Counter)

    def add(
        self,
        observation: Observation | None,
        *,
        dropped: bool = False,
        substitutions: Iterable[tuple[Site, str]] = (),
        flanks: tuple[bool, bool] | None = None,
    ) -> None:
        self.molecules += 1
        if dropped:
            self.dropped += 1
            return
        if observation is None or observation.is_empty:
            self.unusable += 1
            return
        if observation.is_complete:
            self.complete[observation.structure()] += 1
        else:
            self.partial[observation] += 1
        if flanks is not None:
            five, three = flanks
            self.flank_coverage[(observation, five, three)] += 1
            for site, base in substitutions:
                self.substitutions[(observation, site, base)] += 1

    @property
    def has_snp_evidence(self) -> bool:
        return bool(self.flank_coverage)

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
            "substitutions": [
                _observation_dict(o) | {"site": site.key, "base": base, "count": n}
                for (o, site, base), n in self.substitutions.most_common()
            ],
            "flank_coverage": [
                _observation_dict(o) | {"five_prime": five, "three_prime": three, "count": n}
                for (o, five, three), n in self.flank_coverage.most_common()
            ],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> SampleCounts:
        if data.get("schema") not in READABLE:
            raise ValueError(f"unsupported counts schema {data.get('schema')!r}")
        partial: Counter[Observation] = Counter()
        for row in data["partial"]:
            partial[_observation(row)] = row["count"]
        substitutions: Counter[tuple[Observation, Site, str]] = Counter()
        for row in data.get("substitutions", []):
            substitutions[(_observation(row), Site.from_key(row["site"]), row["base"])] = row[
                "count"
            ]
        flank_coverage: Counter[tuple[Observation, bool, bool]] = Counter()
        for row in data.get("flank_coverage", []):
            flank_coverage[(_observation(row), row["five_prime"], row["three_prime"])] = row[
                "count"
            ]
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
            substitutions=substitutions,
            flank_coverage=flank_coverage,
        )

    def write_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2) + "\n")

    @classmethod
    def read_json(cls, path: str | Path) -> SampleCounts:
        return cls.from_dict(json.loads(Path(path).read_text()))


def _observation_dict(o: Observation) -> dict[str, Any]:
    return {
        "observation": o.label,
        "counts": list(o.counts),
        "status": [s.name.lower() for s in o.status],
    }


def _observation(row: dict[str, Any]) -> Observation:
    status = tuple(FieldStatus[s.upper()] for s in row["status"])
    return Observation(tuple(row["counts"]), status)


def _usable(parse: ReadParse) -> Observation | None:
    return parse.observation if parse.usable else None


def _sequence(read: str | FastqRecord) -> tuple[str, str | None]:
    if isinstance(read, FastqRecord):
        return read.sequence, read.quality
    return read, None


def count_reads(
    reads: Iterable[tuple[str | FastqRecord, str | FastqRecord | None]],
    parser: RepeatParser | None = None,
    policy: DiscordancePolicy = DiscordancePolicy.DROP,
    prefer_r1: tuple[bool, ...] = DEFAULT_PREFER_R1,
    *,
    min_base_quality: int = 30,
) -> SampleCounts:
    """Tally ``(r1, r2)`` pairs as sequenced; pass ``r2=None`` for single-end.

    Reads are sequences, or FASTQ records whose base qualities then count towards the SNP
    evidence. A substitution only one read shows needs at least ``min_base_quality``.
    """
    parser = parser or RepeatParser()
    counts = SampleCounts()
    for first_read, second_read in reads:
        r1, q1 = _sequence(first_read)
        first = parser.parse(r1)
        counts.read_outcomes[f"r1.{first.outcome}"] += 1
        first_observation = _usable(first)
        if second_read is None:
            observation = first_observation
            used = [(r1, q1, first)] if first_observation else []
        else:
            sequence, quality = _sequence(second_read)
            r2, q2 = reverse_complement(sequence), quality[::-1] if quality else None
            second = parser.parse(r2)
            counts.read_outcomes[f"r2.{second.outcome}"] += 1
            second_observation = _usable(second)
            joined = join_read_base_pairings(
                first_observation, second_observation, policy=policy, prefer_r1=prefer_r1
            )
            for name, flag in zip(FIELDS, joined.discordant, strict=True):
                if flag:
                    counts.discordant[name] += 1
            if joined.is_discordant and policy is DiscordancePolicy.DROP:
                counts.add(joined.observation, dropped=True)
                continue
            observation = joined.observation
            used = [(r1, q1, first)] if first_observation else []
            if second_observation:
                used.append((r2, q2, second))
        substitutions: frozenset[tuple[Site, str]] = frozenset()
        if observation is not None and any(parse.substitutions for _, _, parse in used):
            substitutions = molecule_substitutions(
                [ReadEvidence(*read) for read in used], observation, min_quality=min_base_quality
            )
        flanks = (
            any(parse.five_anchored for _, _, parse in used),
            any(parse.three_anchored for _, _, parse in used),
        )
        counts.add(observation, substitutions=substitutions, flanks=flanks)
    return counts


def count_fastq(
    r1: str | Path,
    r2: str | Path | None = None,
    parser: RepeatParser | None = None,
    policy: DiscordancePolicy = DiscordancePolicy.DROP,
    prefer_r1: tuple[bool, ...] = DEFAULT_PREFER_R1,
    progress: Progress | None = None,
) -> SampleCounts:
    if r2 is None:
        records = read_fastq(r1, progress)
        pairs: Iterable[tuple[FastqRecord, FastqRecord | None]] = ((rec, None) for rec in records)
    else:
        pairs = read_pairs(r1, r2, progress)
    return count_reads(pairs, parser, policy, prefer_r1)
