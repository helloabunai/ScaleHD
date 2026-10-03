"""Candidate genotypes as mixture models, and their likelihood."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.special import expit, logsumexp

from ...calibration import Stutter, from_transformed, logit
from .kernel import _log_kernel, _log_survival
from .results import Candidate
from .settings import CallerSettings
from .table import _BACKGROUND_SPAN, _Table

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


def _allele_loglik(
    table: _Table, allele: Candidate, stutter: Stutter, params: _Params
) -> np.ndarray:
    """log P(observation | allele) for every row of the table.

    For an allele beyond read length, one row of those per N it could be (L .. L + 59),
    since all its molecules share one N; see _total_loglik.
    """
    s = allele.structure
    ccg_slippage, misread = params.ccg_slippage, params.misread
    exact, lower = table.exact, table.lower

    (exact_values, exact_index), (lower_values, lower_index) = table.cag_exact, table.cag_lower
    unconfirmed_values, unconfirmed_index = table.cag_unconfirmed
    n: np.ndarray | int = s.cag
    if allele.beyond_read_length:
        n = np.arange(s.cag, s.cag + _LONG_ALLELE_SPAN)[:, None]

    def kernel_at(x: np.ndarray) -> np.ndarray:
        return _log_kernel(x, n, stutter, table.window)

    def survival_at(x: np.ndarray) -> np.ndarray:
        return _log_survival(x, n, stutter, table.window)

    # One kernel and one survival call for every value. Each call works out the
    # kernel's normalising sum again, and this runs for every likelihood evaluation.
    split = len(exact_values), len(lower_values)
    kernel, kernel_unconfirmed = np.split(
        kernel_at(np.concatenate([exact_values, unconfirmed_values])), [split[0]], axis=-1
    )
    survival, survival_unconfirmed = np.split(
        survival_at(np.concatenate([lower_values, unconfirmed_values + 1])), [split[1]], axis=-1
    )
    # An unconfirmed CAG count. Likely to be "true", unless a seq error influenced tract end
    # and it goes on.
    error = table.unconfirmed_error
    unconfirmed = np.logaddexp(
        math.log1p(-error) + kernel_unconfirmed, math.log(error) + survival_unconfirmed
    )
    keep, flat = math.log1p(-params.floor), math.log(params.floor)
    span = table.floor_span
    kernel = np.logaddexp(keep + kernel, flat - math.log(span))
    unconfirmed = np.logaddexp(keep + unconfirmed, flat - math.log(span))
    above = np.clip((span - lower_values + 1) / span, 1 / span, 1.0)
    survival = np.logaddexp(keep + survival, flat + np.log(above))
    out = np.zeros((*kernel.shape[:-1], table.size))
    out[..., exact[:, 0]] = kernel[..., exact_index]
    out[..., lower[:, 0]] = survival[..., lower_index]
    out[..., table.unconfirmed_cag] = unconfirmed[..., unconfirmed_index]

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


def _total_loglik(table: _Table, total: np.ndarray) -> float:
    """The sample's log-likelihood from _row_loglik's rows.

    With an allele beyond read length there is a row per N it could be. Its one N is
    averaged over.
    """
    if total.ndim == 1:
        return float(table.fit_weight @ total)
    return float(logsumexp(total @ table.fit_weight)) - math.log(total.shape[0])


def _n_weights(table: _Table, total: np.ndarray) -> np.ndarray | None:
    """How well each N of an allele beyond read length fits, summing to 1 (None if none)."""
    if total.ndim == 1:
        return None
    loglik = total @ table.fit_weight
    return np.asarray(np.exp(loglik - logsumexp(loglik)))


def _per_row(values: np.ndarray, weights: np.ndarray | None) -> np.ndarray:
    """A per-row figure, averaged over N (by weights) if it has a row per N."""
    return values if weights is None or values.ndim == 1 else weights @ values


def _negative_log_posterior(
    theta: np.ndarray, table: _Table, model: _Model, mean: np.ndarray, sd: np.ndarray
) -> float:
    total, _ = _row_loglik(table, model, _unpack(theta, model))
    return -(_total_loglik(table, total) - 0.5 * float((((theta - mean) / sd) ** 2).sum()))
