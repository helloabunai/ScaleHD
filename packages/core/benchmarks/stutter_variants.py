"""Stutter curves unlike the caller's assumptions, for the model benchmarks.

The simulator's default stutter and the caller's stutter priors are the same curve
(``scalehd.calibration.HTT_MISEQ``), so calls on simulated reads partly show the caller
agreeing with the simulator. These curves break that connection; the simulator uses one of
them, the caller keeps its priors.

Only used by ``stutter_robustness.py`` and ``posterior_calibration.py``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scalehd.calibration import HTT_MISEQ, StutterCurve, logit

# (N-1)/N and (N+1)/N stay below this, so an allele is still the tallest peak of its own
# molecules. The caller sizes alleles that way, and so do the truth labels.
_MAX_FIRST = math.log(0.95)
# The step and tail ratios stay within the caller's own bounds.
_RATIO = (logit(1e-3), logit(0.97))
# Columns of a StutterCurve, in order: log(N-1/N), logit(N-2/N-1), logit(tail down),
# log(N+1/N), logit(N+2/N+1), logit(tail up).
_LOG_COLUMNS = (0, 3)
HALF = math.log(0.5)
# Measured spreads of a variety of CAG lengths from scalehd 1.x training data
# need more data hurrhurr
MEASURED_SPREAD = np.array(
    [
        (0.19, 0.23, 0.55, 0.33, 1.14, 1.66),  # CAG 9
        (0.19, 0.23, 0.55, 0.33, 1.14, 1.66),  # CAG 16.7
        (0.23, 0.28, 0.40, 0.31, 1.08, 1.71),  # CAG 22.7
        (0.08, 0.21, 0.27, 0.38, 0.60, 0.39),  # CAG 36.9
        (0.10, 0.19, 0.29, 0.46, 0.55, 0.55),  # CAG 44.5
        (0.10, 0.37, 0.30, 0.14, 0.41, 0.17),  # CAG 51.7
        (0.25, 0.34, 0.10, 0.25, 0.54, 0.42),  # CAG 64.7
    ]
)


def shifted(
    offsets: tuple[float, ...] | np.ndarray = (0.0,) * 6,
    cag_shift: float = 0.0,
    curve: StutterCurve = HTT_MISEQ,
) -> StutterCurve:
    """``curve`` with ``offsets`` added to its six columns (on their log or logit scales),
    and stutter at CAG n behaving as the curve's does at n + ``cag_shift``. An offset can
    be one number per column, or one per point of the curve."""
    columns = (
        curve.log_contraction,
        curve.logit_contraction_step,
        curve.logit_contraction_tail,
        curve.log_expansion,
        curve.logit_expansion_step,
        curve.logit_expansion_tail,
    )
    new = []
    for k, (column, offset) in enumerate(zip(columns, offsets, strict=True)):
        values = np.asarray(column) + offset
        # Log columns are only capped, so N stays the tallest peak; logit ones are bounded.
        top = _MAX_FIRST if k in _LOG_COLUMNS else _RATIO[1]
        bottom = -np.inf if k in _LOG_COLUMNS else _RATIO[0]
        values = np.clip(values, bottom, top)
        new.append(tuple(float(v) for v in values))
    lengths = tuple(float(n - cag_shift) for n in curve.cag_lengths)
    return StutterCurve(lengths, *new, spread=curve.spread)


@dataclass(frozen=True)
class Variant:
    """How the simulator's stutter differs from the caller's priors."""

    name: str
    offsets: tuple[float, ...] = (0.0,) * 6
    cag_shift: float = 0.0
    # Above 0, every sample draws its own offsets from a normal distribution with this
    # many times the caller's prior SD: 1 is what the priors themselves expect.
    spread: float = 0.0
    # Every sample draws one direction per scale, as far at each length as
    # MEASURED_SPREAD says samples vary there.
    measured: bool = False

    def curve(self, draw: int) -> StutterCurve:
        """The simulator's curve; ``draw`` picks a sample's own offsets for a random one."""
        if not any(self.offsets) and not (self.cag_shift or self.spread or self.measured):
            return HTT_MISEQ
        rng = np.random.default_rng([draw, 1234])
        offsets = np.asarray(self.offsets, dtype=float)
        if self.spread:
            offsets = offsets + rng.normal(0.0, self.spread * np.asarray(HTT_MISEQ.spread))
        if self.measured:
            # One column per scale, one row per point of the curve.
            offsets = offsets[:, None] + rng.standard_normal(6)[:, None] * MEASURED_SPREAD.T
        return shifted(offsets, self.cag_shift)


VARIANTS = (
    Variant("as calibrated"),
    Variant("half the stutter", offsets=(HALF, 0, 0, HALF, 0, 0)),
    Variant("double the stutter", offsets=(-HALF, 0, 0, -HALF, 0, 0)),
    Variant("longer tails", offsets=(0, 0, 1, 0, 0, 1)),
    Variant("shorter tails", offsets=(0, 0, -1, 0, 0, -1)),
    Variant("as if 10 CAG longer", cag_shift=10),
    Variant("as if 10 CAG shorter", cag_shift=-10),
    Variant("random, measured spread", measured=True),
    Variant("random, prior spread", spread=1.0),
    Variant("random, twice prior spread", spread=2.0),
)
