"""Fitting each candidate genotype, and scoring it by a Laplace approximation."""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from functools import partial

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from .likelihood import (
    _BALANCE_BOUNDS,
    _NUISANCE_BOUNDS,
    _STUTTER_PARAMS,
    _bounds,
    _Model,
    _negative_log_posterior,
    _Params,
    _prior,
    _unpack,
)
from .settings import CallerSettings
from .table import _Table


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
