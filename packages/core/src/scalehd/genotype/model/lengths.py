"""Settling an allele's exact CAG, and estimating one beyond read length."""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from .fit import _Fit, _hessian
from .likelihood import (
    _LONG_ALLELE_SPAN,
    _STUTTER_BOUNDS,
    _STUTTER_PARAMS,
    _Model,
    _row_loglik,
    _total_loglik,
    _unpack,
)
from .results import Candidate
from .settings import CallerSettings
from .table import _Table

# Posterior at the longest N tried, relative to its peak, above which the reads set no
# upper limit on an allele beyond read length.
_NO_UPPER_LIMIT = 0.01


def _estimate_long_allele(
    table: _Table, fit: _Fit, allele: Candidate, settings: CallerSettings
) -> tuple[int, int, int] | None:
    """Posterior over N for an allele beyond read length (see _length_posterior).

    None when the posterior has not fallen away by the longest N tried: the reads then
    set no upper limit, and any interval would only reflect how far the search went.
    Every N costs a refit, so every 5th is tried first, then each N only where those
    leave any chance.
    """
    index = fit.model.alleles.index(allele)
    lengths = np.arange(allele.structure.cag, allele.structure.cag + _LONG_ALLELE_SPAN)
    coarse = np.unique(np.append(lengths[::5], lengths[-1]))
    chance = _length_posterior(table, fit, index, coarse, settings)
    if chance[-1] > _NO_UPPER_LIMIT * chance.max():
        return None
    likely = coarse[chance > 1e-6 * chance.max()]
    lengths = lengths[(lengths >= likely.min() - 5) & (lengths <= likely.max() + 5)]
    posterior = _length_posterior(table, fit, index, lengths, settings)
    cumulative = np.cumsum(posterior)
    low = int(lengths[int(np.searchsorted(cumulative, 0.05))])
    high = int(lengths[min(int(np.searchsorted(cumulative, 0.95)), lengths.size - 1)])
    return int(lengths[int(np.argmax(posterior))]), low, high


def _local_n(
    table: _Table, fit: _Fit, index: int, settings: CallerSettings, limit: int | None = None
) -> dict[int, float] | None:
    """Posterior over the exact N of one allele, from its own peak region.

    The full model decides which alleles exist. This decides where each one's peak
    is. Only molecules with the allele's own structure and CAG within ``local_radius``
    of it are used, the same molecules for every N compared, with the other allele,
    floor and background held at the global fit. Distant shoulders and junk then cannot
    shift N, which they did when the whole distribution decided it. Alleles within three
    CAG of another with the same structure are left to the joint fit, since their peaks
    overlap, as well as for alleles beyond read length.

    Where reads end near the peak (more than ``unread_share`` of it not read in full),
    that window cuts across the molecules that place N, and picking from it pulled N
    one low. N is then compared on every molecule, as the full model does. N stays
    below ``limit``, where reads stop, if given.
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
    rest = table.exact[:, 1:].all(axis=1) & (table.value[:, 1:] == np.array(s.counts[1:])).all(
        axis=1
    )
    counted = table.exact[:, 0] | table.unconfirmed_cag
    near = rest & counted & (np.abs(table.value[:, 0] - s.cag) <= settings.local_radius)
    if near.sum() < 3:
        return None
    above = table.value[:, 0] >= s.cag - settings.local_radius
    unread = rest & (table.lower[:, 0] | table.unconfirmed_cag) & above
    if table.weight[unread].sum() > settings.unread_share * table.weight[near].sum():
        local = table
    else:
        local = table.subset(near)

    top = s.cag + settings.local_shift
    if limit is not None:
        top = min(top, limit - 1)
    lengths = np.arange(max(1, s.cag - settings.local_shift), top + 1)
    posterior = _length_posterior(local, fit, index, lengths, settings)
    return dict(zip((int(n) for n in lengths), (float(p) for p in posterior), strict=True))


def _length_posterior(
    table: _Table, fit: _Fit, index: int, lengths: np.ndarray, settings: CallerSettings
) -> np.ndarray:
    """Posterior over an allele's N among ``lengths``, as an exact allele of each.

    Its stutter is refit for every N, everything else held at the fit, and each N is
    scored by a Laplace approximation over that refit. Stutter trades off against N
    (more contraction looks like a longer allele -- from my memory so subject to change),
    so holding it fixed claims N far too precisely.
    """
    model = fit.model
    s = model.alleles[index].structure
    k = _STUTTER_PARAMS
    sd = np.array(settings.stutter.spread)
    lower, upper = np.array(_STUTTER_BOUNDS).T
    scores = []
    for n in lengths:
        alleles = list(model.alleles)
        alleles[index] = Candidate(s.with_counts(cag=int(n)))
        trial = _Model(tuple(alleles))
        mean = settings.stutter.transformed(int(n))

        def objective(x: np.ndarray, trial: _Model = trial, mean: np.ndarray = mean) -> float:
            theta = fit.theta.copy()
            theta[k * index : k * (index + 1)] = x
            total, _ = _row_loglik(table, trial, _unpack(theta, trial))
            return -(_total_loglik(table, total) - 0.5 * float((((x - mean) / sd) ** 2).sum()))

        starts = (np.clip(mean, lower, upper), fit.theta[k * index : k * (index + 1)])
        result = min(
            (minimize(objective, x0, method="L-BFGS-B", bounds=_STUTTER_BOUNDS) for x0 in starts),
            key=lambda r: float(r.fun),
        )
        eigenvalues = np.maximum(
            np.linalg.eigvalsh(_hessian(objective, result.x)), 1 / sd.max() ** 2
        )
        scores.append(
            -float(result.fun) - float(np.log(sd).sum()) - 0.5 * float(np.log(eigenvalues).sum())
        )
    values = np.array(scores)
    return np.asarray(np.exp(values - logsumexp(values)))
