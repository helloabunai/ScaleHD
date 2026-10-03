"""A sample's distinct observations as arrays, for the likelihood."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ...counts import SampleCounts
from ...structure import AlleleStructure, FieldStatus, Observation
from .settings import _WINDOW

_EXACT = int(FieldStatus.EXACT)
_LOWER = int(FieldStatus.LOWER_BOUND)
_UNCONFIRMED = int(FieldStatus.UNCONFIRMED)
# Number of values the background component spreads over, per sub-struct
# (cag, caacag, ccgcca, ccg, cct).
_BACKGROUND_SPAN = np.array([250, 4, 4, 25, 5])


# The floor spreads over CAG 1 .. (longest observed + this margin).
_FLOOR_MARGIN = 10


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
        unconfirmed_error: float = 0.03,
    ) -> None:
        self.window = window
        self.unconfirmed_error = unconfirmed_error
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
        unconfirmed = status == _UNCONFIRMED
        self.complete = (status == _EXACT).all(axis=1)
        self.exact = (status == _EXACT) | (unconfirmed & (np.arange(5) > 0))
        self.lower = status == _LOWER
        self.unconfirmed_cag = unconfirmed[:, 0]
        self.total = float(self.weight.sum())
        scale = 1.0 if effective is None or self.total <= 0 else min(1.0, effective / self.total)
        self.fit_weight = self.weight * scale
        self.fit_total = self.total * scale
        # CAG terms depend only on the value, so they are computed once per distinct value.
        self.cag_exact = np.unique(self.value[self.exact[:, 0], 0], return_inverse=True)
        self.cag_lower = np.unique(self.value[self.lower[:, 0], 0], return_inverse=True)
        self.cag_unconfirmed = np.unique(self.value[self.unconfirmed_cag, 0], return_inverse=True)
        self._terms: dict[AlleleStructure, _AlleleTerms] = {}
        seen = self.value[self.exact[:, 0] | self.lower[:, 0] | self.unconfirmed_cag, 0]
        self.floor_span = max(int(seen.max()) if seen.size else 0, 40) + _FLOOR_MARGIN
        span = _BACKGROUND_SPAN
        tail_share = np.clip((span - self.value) / span, 1 / span, 1.0)
        counted = self.exact | (self.unconfirmed_cag[:, None] & (np.arange(5) == 0))
        self.background = (
            np.where(counted, -np.log(span), 0.0) + np.where(self.lower, np.log(tail_share), 0.0)
        ).sum(axis=1)

    def subset(self, mask: np.ndarray) -> _Table:
        """The same table restricted to some rows, sharing the floor span and window."""
        out = object.__new__(_Table)
        out.window, out.floor_span = self.window, self.floor_span
        out.unconfirmed_error = self.unconfirmed_error
        out.observations = [o for o, keep in zip(self.observations, mask, strict=True) if keep]
        out.size = int(mask.sum())
        out.weight, out.fit_weight = self.weight[mask], self.fit_weight[mask]
        out.value, out.exact, out.lower = self.value[mask], self.exact[mask], self.lower[mask]
        out.unconfirmed_cag = self.unconfirmed_cag[mask]
        out.complete, out.background = self.complete[mask], self.background[mask]
        out.total, out.fit_total = float(out.weight.sum()), float(out.fit_weight.sum())
        out.cag_exact = np.unique(out.value[out.exact[:, 0], 0], return_inverse=True)
        out.cag_lower = np.unique(out.value[out.lower[:, 0], 0], return_inverse=True)
        out.cag_unconfirmed = np.unique(out.value[out.unconfirmed_cag, 0], return_inverse=True)
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
