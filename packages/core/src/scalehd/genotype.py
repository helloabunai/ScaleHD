"""Genotype calling from per-molecule repeat structures.

A candidate genotype is a pair of alleles, possibly the same allele twice. Every
molecule in a sample is scored under a mixture::

    P(molecule) = (1 - β) · [w · M(A₁) + (1 - w) · M(A₂)] + β · background

M(A) is what PCR and sequencing make of allele A:

- the CAG length moves by stutter and somatic expansion, following the two-sided
  geometric kernel of `scalehd.calibration`, whose six ratios are fitted per
  allele around length-dependent priors;
- CCG moves by one unit either way with probability γ;
- each of the CAACAG, CCGCCA and CCT counts is misread with probability μ.

A molecule whose CAG tract ran past both reads contributes P(CAG >= its bound), and
unobserved fields are dropped out. An allele beyond read length has
an unknown N >= L and is averaged over N in [L, L + 60), so it must explain the data as
well as a specific N would. The readable contracted molecules then give an estimate of N.

Whether a molecule's CAG tract is read in full depends on its length, so treating
truncated molecules as plain lower bounds would bias N near the read-length limit.
When a sample has many truncated molecules, every CAG count above a threshold T (a
little below their usual bound) is coarsened to "more than T", readable or not, which
makes truncation independent of length.

Each candidate is fitted by maximum a posteriori. Candidates are then compared by a
Laplace approximation to their marginal likelihood, which gives a posterior probability
for the call. A homozygous candidate can only explain a real second peak with stutter
ratios its priors make implausible.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Callable
from dataclasses import dataclass, field
from enum import StrEnum
from functools import partial
from itertools import combinations_with_replacement
from typing import Any

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import expit, logsumexp

from .calibration import HTT_MISEQ, Stutter, StutterCurve, from_transformed, logit
from .counts import SampleCounts
from .structure import AlleleStructure, FieldStatus, Observation

SCHEMA = "scalehd.call/1"

_EXACT = int(FieldStatus.EXACT)
_LOWER = int(FieldStatus.LOWER_BOUND)
# Number of values the background component spreads over, per sub-struct
# (cag, caacag, ccgcca, ccg, cct).
_BACKGROUND_SPAN = np.array([250, 4, 4, 25, 5])
# Optimiser bounds on the transformed scales.
# N-1/N and N+1/N are capped at 1, so an allele is the modal peak of its own molecules,
# as in the usual sizing convention. Without the cap, N can drift one step off the mode
# whenever that lets the tails fit a little better.
_STUTTER_BOUNDS = (
    (math.log(1e-4), 0.0),
    (logit(1e-3), logit(0.97)),
    (logit(1e-3), logit(0.97)),
    (math.log(1e-5), 0.0),
    (logit(1e-3), logit(0.97)),
    (logit(1e-3), logit(0.97)),
)
_STUTTER_PARAMS = len(_STUTTER_BOUNDS)
_BALANCE_BOUNDS = (-6.0, 6.0)
# Range of N averaged over for an allele beyond read length.
_LONG_ALLELE_SPAN = 60
_NUISANCE_BOUNDS = (
    (logit(1e-6), logit(0.2)),  # CCG slippage
    (logit(1e-6), logit(0.2)),  # CCG misassigned
    (logit(1e-6), logit(0.25)),  # misread
    (logit(1e-6), logit(0.3)),  # floor
    (logit(1e-7), logit(0.5)),  # background
)
# The floor spreads over CAG 1 .. (longest observed + this margin).
_FLOOR_MARGIN = 10
# Default reach of the stutter kernel (below, above); see CallerSettings.stutter_window.
_WINDOW = (8, 30)


@dataclass(frozen=True, slots=True)
class CallerSettings:
    stutter: StutterCurve = HTT_MISEQ
    # Reads are PCR copies of a limited number of input templates, not independent
    # molecules, so the likelihood is weighted down to at most this many. Without it,
    # tiny misfits in peak shape outweigh every prior once a sample has 10^5 reads.
    effective_molecules: int | None = 3000
    # How far the stutter kernel reaches (below, above) an allele The peak region is 
    # where the information about N is. With a wide kernel, distant shoulders and junk
    # (common in ScaleHD 1.x alignments) pulled on the tail ratios and tipped N by one 
    # against a clear peak/mode.
    stutter_window: tuple[int, int] = _WINDOW
    # Candidate alleles: the most frequent complete structures, plus CAG +/-1 of the top
    # few.
    max_candidates: int = 6
    neighbour_candidates: int = 3
    min_truncated_fraction: float = 0.03
    # CAG counts above (most common truncation boundary) with some margin.
    truncation_margin: int = 6
    # Candidate genotypes that get a full fit after a quick screen.
    refine: int = 6
    # Exact N of each separated allele is chosen among N +/- local_shift using only its own
    # structure's molecules within +/- local_radius of the peak (re: _local_n).
    local_radius: int = 6
    local_shift: int = 2
    # Priors as (median, SD) on the logit scale. Balance is the shorter allele's share
    # of molecules, whose median in the ScaleHD 1.x training matrix is 0.48.
    balance_prior: tuple[float, float] = (0.0, 0.6)
    ccg_prior: tuple[float, float] = (logit(0.01), 1.0)
    # CCG read as any value at all. In ScaleHD 1.x alignments a few percent of reads that
    # did not span the CCG tract landed on the wrong CCG entirely. Without this term
    # they distort the other allele's stutter fit.
    ccg_misassigned_prior: tuple[float, float] = (logit(0.001), 1.5)
    misread_prior: tuple[float, float] = (logit(0.002), 1.0)
    # Share of an allele's molecules spread flat over every CAG length in its own
    # structure. Without it a long flat tail distorts the stutter fit and shifts peak/N.
    floor_prior: tuple[float, float] = (logit(0.005), 1.5)
    background_prior: tuple[float, float] = (logit(0.002), 1.5)
    # Flag thresholds.
    min_molecules: int = 500
    min_posterior: float = 0.99
    max_background: float = 0.05
    max_dropped: float = 0.15
    balance_range: tuple[float, float] = (0.2, 0.8)
    unexplained_fraction: float = 0.02


class Flag(StrEnum):
    LOW_DEPTH = "low_depth"
    LOW_CONFIDENCE = "low_confidence"
    HOMOZYGOUS = "homozygous"
    # Alleles one CAG apart and otherwise identical: the hardest case to separate from stutter.
    NEIGHBOURING = "neighbouring"
    ATYPICAL = "atypical"
    # No read spanned an allele's CAG tract, so only a lower bound is known.
    BEYOND_READ_LENGTH = "beyond_read_length"
    ALLELE_IMBALANCE = "allele_imbalance"
    HIGH_BACKGROUND = "high_background"
    # A peak the called genotype does not explain: a third allele, contamination or mosaicism.
    UNEXPLAINED_PEAK = "unexplained_peak"
    # Many molecules dropped because their mates disagreed.
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
    # curve beyond the longest CAG it was measured at, so treat it as rough.
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
    # CAG counts above this were treated as "more than" (docstring).
    coarsened_above: int | None = None

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
            "coarsened_above": self.coarsened_above,
        }


@dataclass(frozen=True, slots=True, eq=False)
class _AlleleTerms:
    """Per-row counts that turn the CCG and tail likelihoods into a few multiply-adds."""

    ccg_kept: np.ndarray  # exact CCG equal to the allele's
    ccg_near: np.ndarray  # exact CCG one unit away
    ccg_far: np.ndarray  # exact CCG two or more units away
    ccg_above: np.ndarray  # units by which a CCG lower bound exceeds the allele's
    tail_kept: np.ndarray  # times log(1 - μ)
    tail_swapped: np.ndarray  # times log(μ / 3)


class _Table:
    """The distinct observations of a sample as arrays."""

    def __init__(
        self,
        counts: SampleCounts,
        effective: int | None = None,
        window: tuple[int, int] = _WINDOW,
    ) -> None:
        self.window = window
        rows: list[tuple[Observation, int]] = [
            (Observation.exact(s), n) for s, n in counts.complete.items()
        ]
        rows.extend(counts.partial.items())
        self.observations = [o for o, _ in rows]
        self.size = len(rows)
        self.weight = np.array([n for _, n in rows], dtype=float)
        self.value = np.array([o.counts for o in self.observations], dtype=np.int64).reshape(-1, 5)
        status = np.array(
            [[int(s) for s in o.status] for o in self.observations], dtype=np.int64
        ).reshape(-1, 5)
        self.exact = status == _EXACT
        self.lower = status == _LOWER
        self.complete = self.exact.all(axis=1)
        self.total = float(self.weight.sum())
        scale = 1.0 if effective is None or self.total <= 0 else min(1.0, effective / self.total)
        self.fit_weight = self.weight * scale
        self.fit_total = self.total * scale
        # CAG terms depend only on the value, so they are computed once per distinct value.
        self.cag_exact = np.unique(self.value[self.exact[:, 0], 0], return_inverse=True)
        self.cag_lower = np.unique(self.value[self.lower[:, 0], 0], return_inverse=True)
        self._terms: dict[AlleleStructure, _AlleleTerms] = {}
        seen = self.value[self.exact[:, 0] | self.lower[:, 0], 0]
        self.floor_span = max(int(seen.max()) if seen.size else 0, 40) + _FLOOR_MARGIN
        span = _BACKGROUND_SPAN
        tail_share = np.clip((span - self.value) / span, 1 / span, 1.0)
        self.background = (
            np.where(self.exact, -np.log(span), 0.0) + np.where(self.lower, np.log(tail_share), 0.0)
        ).sum(axis=1)

    def subset(self, mask: np.ndarray) -> _Table:
        """The same table restricted to some rows, sharing the floor span and window."""
        out = object.__new__(_Table)
        out.window, out.floor_span = self.window, self.floor_span
        out.observations = [o for o, keep in zip(self.observations, mask, strict=True) if keep]
        out.size = int(mask.sum())
        out.weight, out.fit_weight = self.weight[mask], self.fit_weight[mask]
        out.value, out.exact, out.lower = self.value[mask], self.exact[mask], self.lower[mask]
        out.complete, out.background = self.complete[mask], self.background[mask]
        out.total, out.fit_total = float(out.weight.sum()), float(out.fit_weight.sum())
        out.cag_exact = np.unique(out.value[out.exact[:, 0], 0], return_inverse=True)
        out.cag_lower = np.unique(out.value[out.lower[:, 0], 0], return_inverse=True)
        out._terms = {}
        return out

    def terms(self, s: AlleleStructure) -> _AlleleTerms:
        if (cached := self._terms.get(s)) is not None:
            return cached
        value, exact, lower = self.value, self.exact, self.lower
        distance = np.abs(value[:, 3] - s.ccg)
        above = lower[:, 3] & (value[:, 3] > s.ccg)
        ccg_kept = (exact[:, 3] & (distance == 0)).astype(float)
        ccg_near = (exact[:, 3] & (distance == 1)).astype(float)
        ccg_far = (exact[:, 3] & (distance >= 2)).astype(float)
        ccg_above = np.where(above, value[:, 3] - s.ccg, 0).astype(float)
        tail_kept = np.zeros(self.size)
        tail_swapped = np.zeros(self.size)
        for k, truth in ((1, s.caacag), (2, s.ccgcca), (4, s.cct)):
            tail_kept += exact[:, k] & (value[:, k] == truth)
            tail_swapped += (exact[:, k] & (value[:, k] != truth)) + (
                lower[:, k] & (value[:, k] > truth)
            )
        terms = _AlleleTerms(ccg_kept, ccg_near, ccg_far, ccg_above, tail_kept, tail_swapped)
        self._terms[s] = terms
        return terms


def _log_heights(delta: np.ndarray, s: Stutter, window: tuple[int, int] = _WINDOW) -> np.ndarray:
    """log peak height relative to N for integer shifts ``delta`` (not normalised).

    Zero outside the stutter window, where only the floor and background reach.
    """
    k = np.abs(delta)
    log_c, log_e = math.log(s.contraction), math.log(s.expansion)
    down = np.where(
        k == 1,
        log_c,
        log_c + math.log(s.contraction_step) + (k - 2) * math.log(s.contraction_tail),
    )
    up = np.where(
        k == 1,
        log_e,
        log_e + math.log(s.expansion_step) + (k - 2) * math.log(s.expansion_tail),
    )
    heights = np.where(delta < 0, down, np.where(delta > 0, up, 0.0))
    inside = (delta >= -window[0]) & (delta <= window[1])
    return np.where(inside, heights, -np.inf)


def _contracted(
    a: np.ndarray, n: np.ndarray, s: Stutter, window: tuple[int, int] = _WINDOW
) -> np.ndarray:
    """Total height of contractions by a .. min(n - 1, window) units (a >= 1)."""
    c, c2, ct = s.contraction, s.contraction_step, s.contraction_tail
    last = np.minimum(n - 1, window[0])
    first = np.where((a <= 1) & (last >= 1), c, 0.0)
    start = np.maximum(a, 2)
    rest = np.where(start <= last, c * c2 * (ct ** (start - 2) - ct ** (last - 1)) / (1 - ct), 0.0)
    return first + rest


def _expanded(b: np.ndarray, s: Stutter, window: tuple[int, int] = _WINDOW) -> np.ndarray:
    """Total height of expansions by b .. window units (b >= 1)."""
    e, e2, et = s.expansion, s.expansion_step, s.expansion_tail
    last = window[1]
    start = np.maximum(b, 2)
    rest = np.where(start <= last, e * e2 * (et ** (start - 2) - et ** (last - 1)) / (1 - et), 0.0)
    return np.where(b <= 1, e, 0.0) + rest


def _log_kernel(
    x: np.ndarray, n: np.ndarray | int, s: Stutter, window: tuple[int, int] = _WINDOW
) -> np.ndarray:
    """log P(CAG = x) for molecules from alleles of n units. x and n broadcast."""
    length = np.asarray(n, dtype=float)
    ones = np.ones_like(length)
    log_z = np.log(1 + _contracted(ones, length, s, window) + _expanded(np.ones(1), s, window))
    return np.where(x >= 1, _log_heights(x - length, s, window) - log_z, -np.inf)


def _log_survival(
    bound: np.ndarray, n: np.ndarray | int, s: Stutter, window: tuple[int, int] = _WINDOW
) -> np.ndarray:
    """log P(CAG ≥ bound) for molecules from alleles of n units. bound and n broadcast."""
    length = np.asarray(n, dtype=float)
    z = 1 + _contracted(np.ones_like(length), length, s, window) + _expanded(np.ones(1), s, window)
    x = np.maximum(bound, 1).astype(float)
    low = x <= length
    below = _contracted(np.where(low, length - x + 1, 1.0), length, s, window) / z
    above = _expanded(np.where(low, 1.0, x - length), s, window) / z
    return np.where(
        low, np.log(np.clip(1 - below, 1e-300, 1.0)), np.log(np.clip(above, 1e-300, None))
    )


def _allele_loglik(
    table: _Table, allele: Candidate, stutter: Stutter, params: _Params
) -> np.ndarray:
    """log P(observation | allele) for every row of the table."""
    s = allele.structure
    ccg_slippage, misread = params.ccg_slippage, params.misread
    exact, lower = table.exact, table.lower
    out = np.zeros(table.size)

    (exact_values, exact_index), (lower_values, lower_index) = table.cag_exact, table.cag_lower
    if allele.beyond_read_length:
        lengths = np.arange(s.cag, s.cag + _LONG_ALLELE_SPAN)[:, None]
        mean = math.log(_LONG_ALLELE_SPAN)
        window = table.window
        kernel = logsumexp(_log_kernel(exact_values[None, :], lengths, stutter, window), axis=0)
        survival = logsumexp(_log_survival(lower_values[None, :], lengths, stutter, window), axis=0)
        kernel, survival = kernel - mean, survival - mean
    else:
        kernel = _log_kernel(exact_values, s.cag, stutter, table.window)
        survival = _log_survival(lower_values, s.cag, stutter, table.window)
    keep, flat = math.log1p(-params.floor), math.log(params.floor)
    span = table.floor_span
    kernel = np.logaddexp(keep + kernel, flat - math.log(span))
    above = np.clip((span - lower_values + 1) / span, 1 / span, 1.0)
    survival = np.logaddexp(keep + survival, flat + np.log(above))
    out[exact[:, 0]] = kernel[exact_index]
    out[lower[:, 0]] = survival[lower_index]

    terms = table.terms(s)
    # CCG: kept, slipped by one unit (γ each way), or misassigned anywhere (ψ).
    psi = params.ccg_misassigned / _BACKGROUND_SPAN[3]
    out += terms.ccg_kept * math.log(1 - 2 * ccg_slippage - params.ccg_misassigned + psi)
    out += terms.ccg_near * math.log(ccg_slippage + psi) + terms.ccg_far * math.log(psi)
    out += terms.ccg_above * math.log(ccg_slippage)
    out += terms.tail_kept * math.log1p(-misread) + terms.tail_swapped * math.log(misread / 3)
    return out


@dataclass(frozen=True, slots=True)
class _Model:
    alleles: tuple[Candidate, ...]  # one if homozygous, else two with the shorter first

    @classmethod
    def of(cls, a: Candidate, b: Candidate) -> _Model:
        if a == b:
            return cls((a,))
        return cls(tuple(sorted((a, b), key=lambda c: c.length_key)))

    @property
    def heterozygous(self) -> bool:
        return len(self.alleles) == 2

    @property
    def label(self) -> str:
        return f"{self.alleles[0].label}/{self.alleles[-1].label}"


@dataclass(frozen=True, slots=True)
class _Params:
    stutter: tuple[Stutter, ...]
    balance: float  # share of the first allele. 1 when homozygous
    ccg_slippage: float
    ccg_misassigned: float
    misread: float
    floor: float
    background: float


def _unpack(theta: np.ndarray, model: _Model) -> _Params:
    n = len(model.alleles)
    k = _STUTTER_PARAMS
    stutter = tuple(from_transformed(theta[k * i : k * (i + 1)]) for i in range(n))
    j = k * n
    balance = 1.0
    if model.heterozygous:
        balance = float(expit(theta[j]))
        j += 1
    gamma, psi, mu, phi, beta = (float(expit(v)) for v in theta[j : j + 5])
    return _Params(stutter, balance, gamma, psi, mu, phi, beta)


def _prior(model: _Model, settings: CallerSettings) -> tuple[np.ndarray, np.ndarray]:
    means: list[float] = []
    sds: list[float] = []
    for allele in model.alleles:
        means.extend(settings.stutter.transformed(allele.structure.cag))
        sds.extend(settings.stutter.spread)
    nuisance = [
        settings.ccg_prior,
        settings.ccg_misassigned_prior,
        settings.misread_prior,
        settings.floor_prior,
        settings.background_prior,
    ]
    if model.heterozygous:
        nuisance.insert(0, settings.balance_prior)
    for mean, sd in nuisance:
        means.append(mean)
        sds.append(sd)
    return np.array(means), np.array(sds)


def _bounds(model: _Model) -> list[tuple[float, float]]:
    bounds = [b for _ in model.alleles for b in _STUTTER_BOUNDS]
    if model.heterozygous:
        bounds.append(_BALANCE_BOUNDS)
    bounds.extend(_NUISANCE_BOUNDS)
    return bounds


def _row_loglik(
    table: _Table, model: _Model, params: _Params
) -> tuple[np.ndarray, list[np.ndarray]]:
    """log P(row) under the full mixture, plus each allele's own log P(row | allele)."""
    per_allele = [
        _allele_loglik(table, allele, stutter, params)
        for allele, stutter in zip(model.alleles, params.stutter, strict=True)
    ]
    mix = per_allele[0]
    if model.heterozygous:
        w = params.balance
        mix = np.logaddexp(math.log(w) + per_allele[0], math.log1p(-w) + per_allele[1])
    total = np.logaddexp(
        math.log1p(-params.background) + mix, math.log(params.background) + table.background
    )
    return total, per_allele


def _negative_log_posterior(
    theta: np.ndarray, table: _Table, model: _Model, mean: np.ndarray, sd: np.ndarray
) -> float:
    total, _ = _row_loglik(table, model, _unpack(theta, model))
    return -(float(table.fit_weight @ total) - 0.5 * float((((theta - mean) / sd) ** 2).sum()))


@dataclass(frozen=True, slots=True, eq=False)
class _Fit:
    model: _Model
    theta: np.ndarray
    params: _Params
    # Laplace approximation to the log marginal likelihood.
    score: float


def _screen(table: _Table, model: _Model, settings: CallerSettings) -> float:
    """Cheap ranking score"""
    mean, sd = _prior(model, settings)
    theta = mean.copy()
    objective = partial(_negative_log_posterior, table=table, model=model, mean=mean, sd=sd)
    if model.heterozygous:
        index = _STUTTER_PARAMS * len(model.alleles)

        def along_balance(z: float) -> float:
            trial = theta.copy()
            trial[index] = z
            return objective(trial)

        best = minimize_scalar(along_balance, bounds=_BALANCE_BOUNDS, method="bounded")
        theta[index] = float(best.x)
    return -objective(theta) - 0.5 * theta.size * math.log(max(table.fit_total, 2.0))


def _hessian(f: Callable[[np.ndarray], float], x: np.ndarray, step: float = 1e-3) -> np.ndarray:
    k = x.size
    hessian = np.empty((k, k))
    centre = f(x)
    unit = np.eye(k) * step
    for i in range(k):
        hessian[i, i] = (f(x + unit[i]) - 2 * centre + f(x - unit[i])) / step**2
        for j in range(i):
            corner = (
                f(x + unit[i] + unit[j])
                - f(x + unit[i] - unit[j])
                - f(x - unit[i] + unit[j])
                + f(x - unit[i] - unit[j])
            )
            hessian[i, j] = hessian[j, i] = corner / (4 * step**2)
    return hessian


def _warm_start(model: _Model, mean: np.ndarray, warm: _Fit) -> np.ndarray:
    """Prior means, with whatever ``warm`` shares with ``model`` copied across."""
    k = _STUTTER_PARAMS
    theta = mean.copy()
    for i, allele in enumerate(model.alleles):
        if allele in warm.model.alleles:
            j = warm.model.alleles.index(allele)
            theta[k * i : k * (i + 1)] = warm.theta[k * j : k * (j + 1)]
    nuisance = len(_NUISANCE_BOUNDS)
    theta[-nuisance:] = warm.theta[-nuisance:]  # shared by every model
    if model.heterozygous and warm.model.heterozygous:
        theta[k * 2] = warm.theta[k * 2]
    return theta


def _refine(
    table: _Table, model: _Model, settings: CallerSettings, warm: _Fit | None = None
) -> _Fit:
    """MAP fit from the prior means and, if given, from a previous fit. keep the better one.

    The background, floor and tail parameters can settle in different locale for
    different candidates, which would make their comparison down to luck.
    """
    mean, sd = _prior(model, settings)
    bounds = _bounds(model)
    lower, upper = np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds])
    objective = partial(_negative_log_posterior, table=table, model=model, mean=mean, sd=sd)
    starts = [np.clip(mean, lower, upper)]
    if warm is not None:
        starts.append(np.clip(_warm_start(model, mean, warm), lower, upper))
    results = [minimize(objective, x0, method="L-BFGS-B", bounds=bounds) for x0 in starts]
    result = min(results, key=lambda r: float(r.fun))
    theta = np.asarray(result.x)

    # log ∫ L·π dθ ≈ log L(θ̂) + log π(θ̂) + (k/2)·log 2π - ½·log|H|. With a normalised
    # Gaussian prior the 2π terms cancel, leaving the expression below. Curvature is
    # floored at the widest prior's, since the posterior can be no flatter than that.
    eigenvalues = np.linalg.eigvalsh(_hessian(objective, theta))
    eigenvalues = np.maximum(eigenvalues, 1.0 / float(sd.max()) ** 2)
    score = -float(result.fun) - float(np.log(sd).sum()) - 0.5 * float(np.log(eigenvalues).sum())
    return _Fit(model, theta, _unpack(theta, model), score)


def candidate_alleles(
    counts: SampleCounts, settings: CallerSettings | None = None
) -> list[Candidate]:
    settings = settings or CallerSettings()
    found: dict[Candidate, None] = {}
    common = [s for s, _ in counts.complete.most_common(settings.max_candidates)]
    for structure in common:
        found[Candidate(structure)] = None
    for structure in common[: settings.neighbour_candidates]:
        for step in (-1, 1):
            if structure.cag + step >= 1:
                found[Candidate(structure.with_counts(cag=structure.cag + step))] = None
    if (long_allele := _long_candidate(counts, settings)) is not None:
        found[long_allele] = None
    return list(found)


def _truncation_bound(counts: SampleCounts, settings: CallerSettings) -> int | None:
    """Most common CAG lower bound, if enough molecules ran past both reads."""
    bound = Counter[int]()
    for observation, n in counts.partial.items():
        if observation.status[0] is FieldStatus.LOWER_BOUND:
            bound[observation.counts[0]] += n
    total = sum(bound.values())
    if not total or total < max(20, settings.min_truncated_fraction * counts.molecules):
        return None
    return bound.most_common(1)[0][0]


def coarsen(counts: SampleCounts, threshold: int) -> SampleCounts:
    """Replace every CAG count above ``threshold`` by the lower bound ``threshold + 1``."""
    out = SampleCounts(
        read_outcomes=counts.read_outcomes,
        discordant=counts.discordant,
        molecules=counts.molecules,
        unusable=counts.unusable,
        dropped=counts.dropped,
    )
    bound = threshold + 1
    for structure, n in counts.complete.items():
        if structure.cag <= threshold:
            out.complete[structure] += n
        else:
            bounded = (FieldStatus.LOWER_BOUND,) + (FieldStatus.EXACT,) * 4
            out.partial[Observation((bound, *structure.counts[1:]), bounded)] += n
    for observation, n in counts.partial.items():
        if observation.status[0] is FieldStatus.UNOBSERVED or observation.counts[0] <= threshold:
            out.partial[observation] += n
        else:
            values = (bound, *observation.counts[1:])
            status = (FieldStatus.LOWER_BOUND, *observation.status[1:])
            out.partial[Observation(values, status)] += n
    return out


def _long_candidate(counts: SampleCounts, settings: CallerSettings) -> Candidate | None:
    if _truncation_bound(counts, settings) is None:
        return None
    bounded = [(o, n) for o, n in counts.partial.items() if o.status[0] is FieldStatus.LOWER_BOUND]
    bound = Counter[int]()
    for observation, n in bounded:
        bound[observation.counts[0]] += n
    defaults = AlleleStructure(1).counts
    rest = []
    for k in range(1, 5):
        seen = Counter[int]()
        for observation, n in bounded:
            if observation.status[k] is FieldStatus.EXACT:
                seen[observation.counts[k]] += n
        rest.append(seen.most_common(1)[0][0] if seen else defaults[k])
    return Candidate(AlleleStructure(bound.most_common(1)[0][0], *rest), beyond_read_length=True)


def _distribution_metrics(
    histogram: np.ndarray, n: int
) -> tuple[float | None, float | None, float | None, float | None]:
    if n >= histogram.size or histogram[n] <= 0:
        return None, None, None, None
    peak = histogram[n]
    backward = float(histogram[max(n - 2, 0) : n].sum() / peak)
    somatic = float(histogram[n + 1 : n + 11].sum() / peak)
    shifts = np.arange(histogram.size) - n
    mass = histogram.sum()
    expansion = float((np.clip(shifts, 0, None) * histogram).sum() / mass)
    contraction = float((np.clip(shifts, None, 0) * histogram).sum() / mass)
    return backward, somatic, expansion, contraction


def _allele_calls(table: _Table, fit: _Fit) -> tuple[AlleleCall, AlleleCall]:
    model, params = fit.model, fit.params
    total, per_allele = _row_loglik(table, model, params)
    shares = [params.balance, 1 - params.balance] if model.heterozygous else [1.0]
    cag = table.value[:, 0]
    exact_cag = table.exact[:, 0]
    size = int(cag.max()) + 41 if table.size else 1
    calls = []
    for allele, stutter, share, own in zip(
        model.alleles, params.stutter, shares, per_allele, strict=True
    ):
        responsibility = np.exp(math.log1p(-params.background) + math.log(share) + own - total)
        attributed = table.weight * responsibility
        histogram = np.bincount(cag[exact_cag], weights=attributed[exact_cag], minlength=size)
        if allele.beyond_read_length:
            backward = somatic = expansion = contraction = None
            estimate = _estimate_long_allele(table, fit, allele)
        else:
            backward, somatic, expansion, contraction = _distribution_metrics(
                histogram, allele.structure.cag
            )
            estimate = None
        calls.append(
            AlleleCall(
                allele,
                share,
                float(attributed.sum()),
                stutter,
                backward,
                somatic,
                expansion,
                contraction,
                estimate,
            )
        )
    if len(calls) == 1:
        half = calls[0]
        half = AlleleCall(
            half.allele,
            0.5,
            half.molecules / 2,
            half.stutter,
            half.backward_slippage,
            half.somatic_mosaicism,
            half.expansion_index,
            half.contraction_index,
            half.cag_estimate,
        )
        return half, half
    return calls[0], calls[1]


def _estimate_long_allele(table: _Table, fit: _Fit, allele: Candidate) -> tuple[int, int, int]:
    """Posterior over N for an allele beyond read length, other parameters held at fit."""
    lengths = np.arange(allele.structure.cag, allele.structure.cag + _LONG_ALLELE_SPAN)
    loglik = []
    for n in lengths:
        exact = Candidate(allele.structure.with_counts(cag=int(n)))
        alleles = tuple(exact if a == allele else a for a in fit.model.alleles)
        total, _ = _row_loglik(table, _Model(alleles), fit.params)
        loglik.append(float(table.fit_weight @ total))
    posterior = np.exp(np.array(loglik) - logsumexp(loglik))
    cumulative = np.cumsum(posterior)
    low = int(lengths[int(np.searchsorted(cumulative, 0.05))])
    high = int(lengths[min(int(np.searchsorted(cumulative, 0.95)), lengths.size - 1)])
    return int(lengths[int(np.argmax(posterior))]), low, high


def _unexplained(table: _Table, fit: _Fit, settings: CallerSettings) -> tuple[tuple[str, int], ...]:
    total, _ = _row_loglik(table, fit.model, fit.params)
    expected = table.total * np.exp(total)
    called = {a.structure for a in fit.model.alleles}
    found = []
    for i in np.flatnonzero(table.complete):
        observed = table.weight[i]
        if observed < settings.unexplained_fraction * table.total:
            continue
        if observed <= 3 * expected[i] + 10:
            continue
        structure = table.observations[i].structure()
        if structure not in called:
            found.append((structure.label, int(observed)))
    return tuple(sorted(found, key=lambda item: -item[1])[:3])


def _local_n(
    table: _Table, fit: _Fit, index: int, settings: CallerSettings
) -> dict[int, float] | None:
    """Posterior over the exact N of one allele, from its own peak region.

    The full model decides which alleles exist. This decides where each one's peak
    is. Only molecules with the allele's own structure and CAG within ``local_radius``
    of it are used, the same molecules for every N compared, with the other allele,
    floor and background held at the global fit. Distant shoulders and junk then cannot
    shift N, which they did when the whole distribution decided it. Alleles within three
    CAG of another with the same structure are left to the joint fit, since their peaks
    overlap, as well as for alleles beyond read length.
    """
    model = fit.model
    allele = model.alleles[index]
    s = allele.structure
    if allele.beyond_read_length:
        return None
    for other in model.alleles:
        if other is allele or other.beyond_read_length:
            continue
        if other.structure.counts[1:] == s.counts[1:] and abs(other.structure.cag - s.cag) <= 3:
            return None
    same_structure = table.complete & (table.value[:, 1:] == np.array(s.counts[1:])).all(axis=1)
    near = same_structure & (np.abs(table.value[:, 0] - s.cag) <= settings.local_radius)
    if near.sum() < 3:
        return None
    local = table.subset(near)

    k = _STUTTER_PARAMS
    sd = np.array(settings.stutter.spread)
    lower, upper = np.array(_STUTTER_BOUNDS).T
    scores: dict[int, float] = {}
    for n in range(max(1, s.cag - settings.local_shift), s.cag + settings.local_shift + 1):
        alleles = list(model.alleles)
        alleles[index] = Candidate(s.with_counts(cag=n))
        trial = _Model(tuple(alleles))
        mean = settings.stutter.transformed(n)

        def objective(x: np.ndarray, trial: _Model = trial, mean: np.ndarray = mean) -> float:
            theta = fit.theta.copy()
            theta[k * index : k * (index + 1)] = x
            total, _ = _row_loglik(local, trial, _unpack(theta, trial))
            return -(float(local.fit_weight @ total) - 0.5 * float((((x - mean) / sd) ** 2).sum()))

        starts = (np.clip(mean, lower, upper), fit.theta[k * index : k * (index + 1)])
        result = min(
            (minimize(objective, x0, method="L-BFGS-B", bounds=_STUTTER_BOUNDS) for x0 in starts),
            key=lambda r: float(r.fun),
        )
        eigenvalues = np.maximum(
            np.linalg.eigvalsh(_hessian(objective, result.x)), 1 / sd.max() ** 2
        )
        scores[n] = (
            -float(result.fun) - float(np.log(sd).sum()) - 0.5 * float(np.log(eigenvalues).sum())
        )
    values = np.array(list(scores.values()))
    posterior = np.exp(values - logsumexp(values))
    return dict(zip(scores, (float(p) for p in posterior), strict=True))


def _same_configuration(a: _Model, b: _Model, shift: int) -> bool:
    """Same alleles apart from CAG differences of at most ``shift``."""
    if len(a.alleles) != len(b.alleles):
        return False
    for x, y in zip(a.alleles, b.alleles, strict=True):
        if x.beyond_read_length != y.beyond_read_length:
            return False
        if x.structure.counts[1:] != y.structure.counts[1:]:
            return False
        if abs(x.structure.cag - y.structure.cag) > shift:
            return False
    return True


def call_genotype(counts: SampleCounts, settings: CallerSettings | None = None) -> GenotypeCall:
    settings = settings or CallerSettings()
    candidates = candidate_alleles(counts, settings)
    threshold = None
    if (bound := _truncation_bound(counts, settings)) is not None:
        threshold = bound - settings.truncation_margin
        counts = coarsen(counts, threshold)
    table = _Table(counts, settings.effective_molecules, settings.stutter_window)
    if table.total <= 0:
        raise NoMoleculesError("no usable molecules to genotype")

    models = sorted(
        {_Model.of(a, b) for a, b in combinations_with_replacement(candidates, 2)},
        key=lambda m: m.label,
    )
    screened = sorted(models, key=lambda m: -_screen(table, m, settings))
    keep = screened[: settings.refine]
    for wanted in (True, False):  # always refine the best of each zygosity
        best = next((m for m in screened if m.heterozygous is wanted), None)
        if best is not None and best not in keep:
            keep.append(best)

    fits = [_refine(table, model, settings) for model in keep]
    leader = max(fits, key=lambda f: f.score)
    fits = [f if f is leader else _refine(table, f.model, settings, warm=leader) for f in fits]
    scores = np.array([f.score for f in fits])
    log_posterior = scores - logsumexp(scores)
    order = np.argsort(-log_posterior)
    best_fit = fits[int(order[0])]

    # Settle each allele's exact N locally, then refit the refined genotype.
    refined = list(best_fit.model.alleles)
    local_posteriors: list[dict[int, float] | None] = []
    for i, allele in enumerate(best_fit.model.alleles):
        by_n = _local_n(table, best_fit, i, settings)
        local_posteriors.append(by_n)
        if by_n:
            refined[i] = Candidate(
                allele.structure.with_counts(cag=max(by_n, key=by_n.__getitem__))
            )
    configuration = best_fit.model
    if tuple(refined) != best_fit.model.alleles:
        best_fit = _refine(table, _Model(tuple(refined)), settings, warm=best_fit)

    # P(call) = P(this configuration of alleles) x P(each allele's exact N | configuration).
    same = np.array(
        [_same_configuration(f.model, configuration, settings.local_shift) for f in fits]
    )
    log_configuration = float(logsumexp(log_posterior[same]))
    log_n = [
        math.log(max(p[n.structure.cag], 1e-300))
        for p, n in zip(local_posteriors, refined, strict=True)
        if p
    ]
    log_call = log_configuration + sum(log_n)
    posterior = min(1.0, float(np.exp(log_call)))
    alternatives_found: list[tuple[str, float]] = [
        (fits[int(i)].model.label, float(np.exp(log_posterior[i]))) for i in order if not same[i]
    ]
    for i, by_n in enumerate(local_posteriors):
        if not by_n:
            continue
        for n, p in by_n.items():
            if n == refined[i].structure.cag:
                continue
            variant = list(refined)
            variant[i] = Candidate(refined[i].structure.with_counts(cag=n))
            share = float(np.exp(log_call)) / max(by_n[refined[i].structure.cag], 1e-300) * p
            alternatives_found.append((_Model(tuple(variant)).label, share))
    alternatives = tuple(sorted(alternatives_found, key=lambda a: -a[1])[:3])
    error = max(1.0 - posterior, sum(p for _, p in alternatives_found), 1e-10)
    quality = min(99.0, -10 * math.log10(error))

    alleles = _allele_calls(table, best_fit)
    unexplained = _unexplained(table, best_fit, settings)
    params = best_fit.params
    flags = _flags(counts, best_fit, alleles, posterior, unexplained, settings)
    return GenotypeCall(
        alleles=alleles,
        posterior=posterior,
        quality=quality,
        alternatives=alternatives,
        flags=flags,
        molecules=int(table.total),
        background=params.background,
        ccg_slippage=params.ccg_slippage,
        misread=params.misread,
        unexplained=unexplained,
        coarsened_above=threshold,
    )


def _flags(
    counts: SampleCounts,
    fit: _Fit,
    alleles: tuple[AlleleCall, AlleleCall],
    posterior: float,
    unexplained: tuple[tuple[str, int], ...],
    settings: CallerSettings,
) -> tuple[Flag, ...]:
    flags = []
    if counts.molecules - counts.unusable - counts.dropped < settings.min_molecules:
        flags.append(Flag.LOW_DEPTH)
    if posterior < settings.min_posterior:
        flags.append(Flag.LOW_CONFIDENCE)
    first, second = (a.allele for a in alleles)
    if first == second:
        flags.append(Flag.HOMOZYGOUS)
    elif (
        not first.beyond_read_length
        and not second.beyond_read_length
        and abs(first.structure.cag - second.structure.cag) == 1
        and first.structure.counts[1:] == second.structure.counts[1:]
    ):
        flags.append(Flag.NEIGHBOURING)
    if any(not a.allele.structure.is_typical for a in alleles):
        flags.append(Flag.ATYPICAL)
    if any(a.allele.beyond_read_length for a in alleles):
        flags.append(Flag.BEYOND_READ_LENGTH)
    low, high = settings.balance_range
    if fit.model.heterozygous and not low <= fit.params.balance <= high:
        flags.append(Flag.ALLELE_IMBALANCE)
    if fit.params.background > settings.max_background:
        flags.append(Flag.HIGH_BACKGROUND)
    if unexplained:
        flags.append(Flag.UNEXPLAINED_PEAK)
    if counts.molecules and counts.dropped / counts.molecules > settings.max_dropped:
        flags.append(Flag.HIGH_DISCORDANCE)
    return tuple(flags)
