"""Repeat-structure model for the HTT exon 1 CAG/CCG region.

Along the forward (CAG) strand an allele is five tandem tracts, typically in this order::

    (CAG)n (CAACAG)a (CCGCCA)b (CCG)m (CCT)k

and is labelled ``n_a_b_m_k``, matching ScaleHD 1.x reference names. The common
allele is ``n_1_1_7_2``. Known atypical forms fit the same structure. loss of the CAA
interruption is ``n+2_0_1_m_k``, duplication of CAACAG is ``n_2_1_m_k``.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from enum import IntEnum

UNITS: tuple[str, ...] = ("CAG", "CAACAG", "CCGCCA", "CCG", "CCT")
FIELDS: tuple[str, ...] = ("cag", "caacag", "ccgcca", "ccg", "cct")

type Counts = tuple[int, int, int, int, int]


@dataclass(frozen=True, slots=True, order=True)
class AlleleStructure:
    cag: int
    caacag: int = 1
    ccgcca: int = 1
    ccg: int = 7
    cct: int = 2

    def __post_init__(self) -> None:
        for name in FIELDS:
            if getattr(self, name) < 0:
                raise ValueError(f"{name} count must be non-negative")

    @classmethod
    def from_label(cls, label: str) -> AlleleStructure:
        parts = label.split("_")
        if len(parts) != len(FIELDS) or not all(p.isdigit() for p in parts):
            raise ValueError(f"expected a label like '17_1_1_7_2', got {label!r}")
        return cls(*(int(p) for p in parts))

    @classmethod
    def from_counts(cls, counts: Sequence[int]) -> AlleleStructure:
        if len(counts) != len(FIELDS):
            raise ValueError(f"expected {len(FIELDS)} counts, got {len(counts)}")
        return cls(*counts)

    @property
    def counts(self) -> Counts:
        return (self.cag, self.caacag, self.ccgcca, self.ccg, self.cct)

    @property
    def label(self) -> str:
        return "_".join(str(c) for c in self.counts)

    @property
    def is_typical(self) -> bool:
        """Common intervening sequence and CCT tract as per scalehd 1.x."""
        return self.caacag == 1 and self.ccgcca == 1 and self.cct == 2

    @property
    def polyglutamine_length(self) -> int:
        """CAA and CAG codons in the tract."""
        return self.cag + 2 * self.caacag

    def repeat_sequence(self) -> str:
        return "".join(unit * n for unit, n in zip(UNITS, self.counts, strict=True))

    def with_counts(self, **changes: int) -> AlleleStructure:
        values = dict(zip(FIELDS, self.counts, strict=True)) | changes
        return AlleleStructure(**values)

    def __str__(self) -> str:
        return self.label


class FieldStatus(IntEnum):
    UNOBSERVED = 0
    # The read ended inside this tract: it is at least this long.
    LOWER_BOUND = 1
    EXACT = 2
    # Tract ended too close to end of read = caution
    UNCONFIRMED = 3


@dataclass(frozen=True, slots=True)
class Observation:
    """Repeat counts seen in a read or read pair.

    A read that ends inside the repeat only gives a lower bound for the tract it
    ends in, and nothing for the tracts beyond it (very long CAG e.g.). One that ends
    just past a tract's end gives that tract's count, unconfirmed.
    """

    counts: Counts
    status: tuple[FieldStatus, ...]

    @classmethod
    def exact(cls, structure: AlleleStructure) -> Observation:
        return cls(structure.counts, (FieldStatus.EXACT,) * len(FIELDS))

    @property
    def is_complete(self) -> bool:
        return all(s is FieldStatus.EXACT for s in self.status)

    @property
    def is_empty(self) -> bool:
        return all(s is FieldStatus.UNOBSERVED for s in self.status)

    def structure(self) -> AlleleStructure:
        if not self.is_complete:
            raise ValueError(f"observation {self.label} is incomplete")
        return AlleleStructure.from_counts(self.counts)

    @property
    def label(self) -> str:
        """Like ``n_a_b_m_k``, with ``n+`` for lower bounds, ``n~`` for unconfirmed and
        ``?`` for unobserved."""
        parts = []
        for count, status in zip(self.counts, self.status, strict=True):
            match status:
                case FieldStatus.EXACT:
                    parts.append(str(count))
                case FieldStatus.LOWER_BOUND:
                    parts.append(f"{count}+")
                case FieldStatus.UNCONFIRMED:
                    parts.append(f"{count}~")
                case FieldStatus.UNOBSERVED:
                    parts.append("?")
        return "_".join(parts)

    def __str__(self) -> str:
        return self.label
