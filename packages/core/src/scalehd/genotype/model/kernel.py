"""The PCR stutter kernel: how molecules' CAG counts spread around their allele's."""

from __future__ import annotations

import math

import numpy as np

from ...calibration import Stutter
from .settings import _WINDOW


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
    """log P(CAG >= bound) for molecules from alleles of n units. bound and n broadcast."""
    length = np.asarray(n, dtype=float)
    z = 1 + _contracted(np.ones_like(length), length, s, window) + _expanded(np.ones(1), s, window)
    x = np.maximum(bound, 1).astype(float)
    low = x <= length
    below = _contracted(np.where(low, length - x + 1, 1.0), length, s, window) / z
    above = _expanded(np.where(low, 1.0, x - length), s, window) / z
    return np.where(
        low, np.log(np.clip(1 - below, 1e-300, 1.0)), np.log(np.clip(above, 1e-300, None))
    )
