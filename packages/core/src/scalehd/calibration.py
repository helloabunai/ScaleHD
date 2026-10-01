"""Length-dependent PCR stutter, shared by the simulator and the genotype caller.

Stutter around an allele of N CAG units is described by peak-height ratios, three on
each side:

- contraction: (N-1)/N, then (N-2)/(N-1), then one decay for every further step
- expansion: (N+1)/N, then (N+2)/(N+1), then one decay for every further step

The first two steps get their own ratios because real expanded alleles have a sharp
N+1 step followed by a long, flatter somatic tail, which one geometric cannot follow.
In blood DNA the expansion side includes ordinary somatic expansion.

The default curve is the per-length median of these ratios around the labelled alleles
in the ScaleHD 1.x training matrix (``legacy/ScaleHD/train/raw_matrix.csv``, 594 MiSeq
samples aligned by ScaleHD 1.x). Alleles closer than eight CAG to their partner in the
same CCG were dropped, and so were alleles with fewer than 500 reads at N.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np


def logit(p: float) -> float:
    return float(np.log(p) - np.log1p(-p))


def expit(x: float) -> float:
    return float(1.0 / (1.0 + np.exp(-x)))


@dataclass(frozen=True, slots=True)
class Stutter:
    contraction: float  # (N-1)/N
    contraction_step: float  # (N-2)/(N-1)
    contraction_tail: float  # (N-k-1)/(N-k) for k >= 2
    expansion: float  # (N+1)/N
    expansion_step: float  # (N+2)/(N+1)
    expansion_tail: float  # (N+k+1)/(N+k) for k >= 2

    def heights(self, shifts: np.ndarray) -> np.ndarray:
        """Peak heights relative to N for the given CAG shifts (not normalised)."""
        k = np.abs(shifts)
        down = np.where(
            k == 1,
            self.contraction,
            self.contraction
            * self.contraction_step
            * self.contraction_tail ** np.maximum(k - 2, 0),
        )
        up = np.where(
            k == 1,
            self.expansion,
            self.expansion * self.expansion_step * self.expansion_tail ** np.maximum(k - 2, 0),
        )
        return np.where(shifts < 0, down, np.where(shifts > 0, up, 1.0))

    def kernel(self, cag: int, max_shift: int = 40) -> tuple[np.ndarray, np.ndarray]:
        """CAG shifts and their probabilities for a template of ``cag`` units."""
        shifts = np.arange(max(-max_shift, 1 - cag), max_shift + 1)
        weights = self.heights(shifts)
        return shifts, weights / weights.sum()

    def as_dict(self) -> dict[str, float]:
        return {
            "contraction": self.contraction,
            "contraction_step": self.contraction_step,
            "contraction_tail": self.contraction_tail,
            "expansion": self.expansion,
            "expansion_step": self.expansion_step,
            "expansion_tail": self.expansion_tail,
        }


@dataclass(frozen=True, slots=True)
class StutterCurve:
    """Stutter ratios as a function of CAG length, interpolated and made up.

    ``cag_lengths`` are the CAG lengths the ratios were measured at. In between, values
    are interpolated linearly; beyond the shortest and longest they are at a constant.
    """

    cag_lengths: tuple[float, ...]
    log_contraction: tuple[float, ...]
    logit_contraction_step: tuple[float, ...]
    logit_contraction_tail: tuple[float, ...]
    log_expansion: tuple[float, ...]
    logit_expansion_step: tuple[float, ...]
    logit_expansion_tail: tuple[float, ...]
    # Prior widths on the same scales, used by the caller. roughly the between-sample
    # spread in the training matrix, widened to allow for other PCR conditions.
    spread: tuple[float, ...] = (0.35, 0.45, 0.8, 0.6, 1.0, 1.2)

    def transformed(self, cag: float) -> np.ndarray:
        columns = (
            self.log_contraction,
            self.logit_contraction_step,
            self.logit_contraction_tail,
            self.log_expansion,
            self.logit_expansion_step,
            self.logit_expansion_tail,
        )
        return np.array([np.interp(cag, self.cag_lengths, column) for column in columns])

    def at(self, cag: float) -> Stutter:
        return from_transformed(self.transformed(cag))


def from_transformed(values: np.ndarray) -> Stutter:
    c1, c2, ct, e1, e2, et = (float(v) for v in values)
    return Stutter(np.exp(c1), expit(c2), expit(ct), np.exp(e1), expit(e2), expit(et))


# Medians from the ScaleHD 1.x training matrix at the mean CAG of each length bin.
HTT_MISEQ = StutterCurve(
    cag_lengths=(9.0, 16.7, 22.7, 36.9, 44.5, 51.7, 64.7),
    log_contraction=(-2.3, -1.73, -1.32, -0.62, -0.47, -0.35, -0.19),
    logit_contraction_step=(-2.3, -1.86, -1.42, -0.69, -0.38, 0.12, 0.71),
    logit_contraction_tail=(-0.5, -0.32, -0.08, -0.07, -0.03, 0.59, 0.70),
    log_expansion=(-4.6, -4.25, -3.82, -2.33, -1.24, -0.44, -0.15),
    logit_expansion_step=(-2.4, -2.37, -2.33, -1.40, -1.15, 0.02, 1.08),
    logit_expansion_tail=(0.0, 0.0, 0.0, 0.40, -0.18, 0.13, 0.64),
)


@lru_cache(maxsize=4096)
def cag_kernel(curve: StutterCurve, cag: int) -> tuple[np.ndarray, np.ndarray]:
    return curve.at(cag).kernel(cag)
