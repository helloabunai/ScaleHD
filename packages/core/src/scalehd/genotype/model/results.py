"""What the caller reports: the call, each allele's figures, and flags."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any

from ...calibration import Stutter
from ...structure import AlleleStructure

SCHEMA = "scalehd.call/2"


class Flag(StrEnum):
    LOW_DEPTH = "low_depth"
    LOW_CONFIDENCE = "low_confidence"
    HOMOZYGOUS = "homozygous"
    # Alleles one CAG apart and otherwise identical: the hardest case to separate from stutter.
    NEIGHBOURING = "neighbouring"
    # Alleles further apart, but close enough that one's stutter is a fair share of the
    # other's peak, so their exact CAG are hard to tell apart. Most often two long alleles.
    # e.g. something like 35 / 37
    CLOSE_ALLELES = "close_alleles"
    ATYPICAL = "atypical"
    # No read spanned an allele's CAG tract, so only a lower bound is known.
    BEYOND_READ_LENGTH = "beyond_read_length"
    ALLELE_IMBALANCE = "allele_imbalance"
    HIGH_BACKGROUND = "high_background"
    # A peak the called genotype does not explain: a third allele, contamination or mosaicism.
    UNEXPLAINED_PEAK = "unexplained_peak"
    # Many molecules dropped because their read base pairings disagreed.
    HIGH_DISCORDANCE = "high_discordance"


class NoMoleculesError(ValueError):
    pass


@dataclass(frozen=True, slots=True, order=True)
class Candidate:
    """A possible allele.

    ``beyond_read_length`` means no read spanned its CAG tract, so its CAG count is at
    least ``structure.cag``.
    """

    structure: AlleleStructure
    beyond_read_length: bool = False

    @property
    def label(self) -> str:
        if not self.beyond_read_length:
            return self.structure.label
        _, rest = self.structure.label.split("_", 1)
        return f"{self.structure.cag}+_{rest}"

    @property
    def length_key(self) -> tuple[int, int]:
        """Orders alleles by amplicon length; one beyond read length sorts last."""
        return (int(self.beyond_read_length), len(self.structure.repeat_sequence()))


@dataclass(frozen=True, slots=True)
class AlleleCall:
    allele: Candidate
    fraction: float
    molecules: float
    stutter: Stutter
    # ScaleHD 1.x metrics: (N-2 + N-1)/N and (N+1 .. N+10)/N.
    backward_slippage: float | None
    somatic_mosaicism: float | None
    # Mean shift of the allele's molecules from N, split into the expanded and contracted
    # sides (the contraction index is negative). PCR stutter contributes to both.
    expansion_index: float | None
    contraction_index: float | None
    # For an allele beyond read length, a model-based estimate of N from its readable contractions
    # (most likely, lower, upper) bounds of a 90% interval. It leans on the stutter
    # curve beyond the longest CAG it was measured at, so treat it as rough. None when
    # the reads set no upper limit.
    cag_estimate: tuple[int, int, int] | None = None

    @property
    def label(self) -> str:
        return self.allele.label

    def to_dict(self) -> dict[str, Any]:
        s = self.allele.structure
        return {
            "structure": self.label,
            "beyond_read_length": self.allele.beyond_read_length,
            "cag": s.cag,
            "caacag": s.caacag,
            "ccgcca": s.ccgcca,
            "ccg": s.ccg,
            "cct": s.cct,
            "typical": s.is_typical,
            "polyglutamine_length": None
            if self.allele.beyond_read_length
            else s.polyglutamine_length,
            "fraction": self.fraction,
            "molecules": self.molecules,
            "stutter": self.stutter.as_dict(),
            "backward_slippage": self.backward_slippage,
            "somatic_mosaicism": self.somatic_mosaicism,
            "expansion_index": self.expansion_index,
            "contraction_index": self.contraction_index,
            "cag_estimate": list(self.cag_estimate) if self.cag_estimate else None,
        }


@dataclass(frozen=True, slots=True)
class GenotypeCall:
    alleles: tuple[AlleleCall, AlleleCall]
    posterior: float
    # Phred-scaled probability that the call is wrong, capped at 99.
    quality: float
    alternatives: tuple[tuple[str, float], ...]
    flags: tuple[Flag, ...]
    molecules: int
    background: float
    ccg_slippage: float
    misread: float
    unexplained: tuple[tuple[str, int], ...] = field(default=())

    @property
    def label(self) -> str:
        return "/".join(a.label for a in self.alleles)

    @property
    def homozygous(self) -> bool:
        return self.alleles[0].allele == self.alleles[1].allele

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "genotype": self.label,
            "posterior": self.posterior,
            "quality": self.quality,
            "flags": [str(f) for f in self.flags],
            "alleles": [a.to_dict() for a in self.alleles],
            "alternatives": [{"genotype": g, "posterior": p} for g, p in self.alternatives],
            "molecules": self.molecules,
            "background": self.background,
            "ccg_slippage": self.ccg_slippage,
            "misread": self.misread,
            "unexplained": [{"structure": s, "molecules": n} for s, n in self.unexplained],
        }
