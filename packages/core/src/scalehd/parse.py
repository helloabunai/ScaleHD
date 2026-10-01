"""Read-level parsing of the HTT repeat region.

A read is placed by short anchors at the flank ends next to the repeat. The bases
between the anchors are split into established tract structure of `scalehd.structure`
by a small dynamic script that tolerates substitutions and single-base ins/dels.

Reads must be in forward (CAG-strand) orientation, so reverse-complement R2 first.
A read that ends inside the repeat has one "open" end, where a partial unit or a
little junk is allowed and the tract it ends in only gives lower trust.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, replace
from enum import StrEnum
from functools import lru_cache
from typing import Literal

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from .amplicon import HTT_AMPLICON, AmpliconSpec
from .structure import UNITS, Counts, FieldStatus, Observation


class ReadOutcome(StrEnum):
    SPANNING = "spanning"  # both anchors found; the whole repeat is in the read
    TRUNCATED = "truncated"  # one anchor found; the read ends inside the repeat
    NO_ANCHOR = "no_anchor"  # off-target, primer dimer or junk?
    NONCONFORMING = "nonconforming"  # too many differences from the established structure


@dataclass(frozen=True, slots=True)
class ParserSettings:
    anchor_length: int = 20
    anchor_max_mismatches: int = 2
    # Shortest anchor fragment accepted where a read ends inside an anchor; must match exactly.
    min_partial_anchor: int = 10
    substitution_cost: int = 2
    indel_cost: int = 3
    # Per-base cost of soft-clipping at an open end. at most anchor_length bases are left over.
    clip_cost: int = 1
    # At an open end a tract boundary is only trusted when at least this many error-free
    # bases follow it. One substitution in a read's tail
    # can fake a boundary (CAG>CAA reads as CAACAG) but rarely gives a run this long.
    confirm_bases: int = 12
    # A read whose decomposition needs more than max_errors + error_rate * len(region)
    # substitution-equivalents is reported as nonconforming.
    max_errors: float = 2.0
    error_rate: float = 0.03


@dataclass(frozen=True, slots=True)
class ReadParse:
    outcome: ReadOutcome
    observation: Observation | None = None
    errors: float = 0.0
    region: str = ""

    @property
    def usable(self) -> bool:
        return (
            self.outcome in (ReadOutcome.SPANNING, ReadOutcome.TRUNCATED)
            and self.observation is not None
            and not self.observation.is_empty
        )


_EXACT = re.compile("".join(f"((?:{unit})*)" for unit in UNITS))
_INF = 1 << 30


def find_anchor(
    read: str,
    anchor: str,
    *,
    max_mismatches: int,
    start: int = 0,
    min_partial: int | None = None,
    partial_side: Literal["left", "right"] | None = None,
) -> int | None:
    """Return the start of ``anchor`` in ``read[start:]``, in read coordinates.

    With ``partial_side`` the anchor may overhang that end of the read, as long as at
    least ``min_partial`` bases overlap and match exactly. The returned position can
    then be negative ("left") or leave the anchor running past the end ("right").
    """
    hit = read.find(anchor, start)
    if hit >= 0:
        return hit

    width = len(anchor)
    overhang = width - min_partial if (partial_side and min_partial) else 0
    left = overhang if partial_side == "left" else 0
    right = overhang if partial_side == "right" else 0
    body = np.frombuffer(read[start:].encode("ascii"), dtype=np.uint8)
    padded = np.concatenate((np.zeros(left, np.uint8), body, np.zeros(right, np.uint8)))
    if padded.size < width:
        return None

    windows = sliding_window_view(padded, width)
    present = windows != 0
    target = np.frombuffer(anchor.encode("ascii"), dtype=np.uint8)
    mismatches = ((windows != target) & present).sum(axis=1)
    overlap = present.sum(axis=1)
    full = overlap == width
    ok = (overlap >= (min_partial or width)) & (mismatches <= np.where(full, max_mismatches, 0))
    if not ok.any():
        return None
    # Prefer full-length hits, then fewer mismatches; ties go to the left-most hit.
    score = np.where(full, mismatches, max_mismatches + 1)
    candidates = np.flatnonzero(ok)
    best = int(candidates[np.argmin(score[candidates])])
    return start + best - left


@lru_cache(maxsize=1 << 16)
def _unit_cost(chunk: str, unit: str, substitution: int, indel: int) -> int:
    """Weighted edit distance between a read chunk and one repeat unit."""
    if chunk == unit:
        return 0
    previous = [j * indel for j in range(len(unit) + 1)]
    for i, base in enumerate(chunk, 1):
        current = [i * indel]
        for j, expected in enumerate(unit, 1):
            current.append(
                min(
                    previous[j] + indel,
                    current[j - 1] + indel,
                    previous[j - 1] + (0 if base == expected else substitution),
                )
            )
        previous = current
    return previous[-1]


@dataclass(frozen=True, slots=True)
class Decomposition:
    counts: Counts
    # Substitution-equivalents, excluding soft-clipped bases.
    errors: float
    # Every whole unit in read order, as (tract index, bases consumed, cost).
    units: tuple[tuple[int, int, int], ...]
    # Bases left over at an open start or end: a partial unit or soft-clipped junk.
    head: str = ""
    tail: str = ""
    # Fields whose count changes depending on how equal-cost ties are broken, as when
    # a substitution leaves one unit equally close to CCG and CCT.
    ambiguous: tuple[bool, ...] = (False,) * len(UNITS)

    @property
    def path(self) -> tuple[int, ...]:
        return tuple(state for state, _, _ in self.units)


def _edge_cost(piece: str, states: range, *, at_end: bool, sub: int, clip: int) -> tuple[int, bool]:
    """Cost of leftover bases at an open end, and whether they were soft-clipped.

    The leftover is either part of a unit from ``states`` (a prefix of it at the
    end of a read, a suffix at the start), costed by mismatches, or clipped.
    """
    if not piece:
        return 0, False
    best, clipped = clip * len(piece), True
    for k in states:
        unit = UNITS[k]
        if len(piece) < len(unit):
            expected = unit[: len(piece)] if at_end else unit[-len(piece) :]
            total = sub * sum(a != b for a, b in zip(piece, expected, strict=True))
            if total < best:
                best, clipped = total, False
    return best, clipped


@lru_cache(maxsize=1 << 18)
def decompose(
    region: str, open_start: bool, open_end: bool, settings: ParserSettings
) -> Decomposition | None:
    """Split ``region`` into whole tract units, allowing one open (truncated) end."""
    if open_start and open_end:
        raise ValueError("a region cannot be open at both ends")
    if not (open_start or open_end) and (match := _EXACT.fullmatch(region)):
        found = tuple(len(g) // len(u) for g, u in zip(match.groups(), UNITS, strict=True))
        whole = tuple((k, len(UNITS[k]), 0) for k in range(len(UNITS)) for _ in range(found[k]))
        return Decomposition(found, 0.0, whole)  # type: ignore[arg-type]

    early = _dynamic_parse(region, open_start, open_end, settings, late=False)
    if early is None or not _error_at_boundary(early):
        return early
    late = _dynamic_parse(region, open_start, open_end, settings, late=True)
    if late is None or late.counts == early.counts:
        return early
    ambiguous = tuple(a != b for a, b in zip(early.counts, late.counts, strict=True))
    return replace(early, ambiguous=ambiguous)


def _error_at_boundary(d: Decomposition) -> bool:
    """Whether a unit with errors sits next to a change of tract."""
    units = d.units
    for index, (state, _, unit_cost) in enumerate(units):
        if not unit_cost:
            continue
        if index > 0 and units[index - 1][0] != state:
            return True
        if index + 1 < len(units) and units[index + 1][0] != state:
            return True
    return False


def _dynamic_parse(
    region: str, open_start: bool, open_end: bool, settings: ParserSettings, *, late: bool
) -> Decomposition | None:
    """Cheapest decomposition. on ties, tract changes go as early or as late as possible."""
    states = len(UNITS)
    length = len(region)
    sub, indel, clip = settings.substitution_cost, settings.indel_cost, settings.clip_cost
    max_edge = min(settings.anchor_length, length)

    # cost[i][k]: cheapest parse of region[:i] that has reached tract k.
    # back[i][k]: (previous i, previous k, whether a unit was consumed); None marks a start.
    cost = [[_INF] * states for _ in range(length + 1)]
    back: list[list[tuple[int, int, bool] | None]] = [[None] * states for _ in range(length + 1)]
    if open_start:
        for c in range(max_edge + 1):
            for k in range(states):
                cost[c][k] = _edge_cost(region[:c], range(k + 1), at_end=False, sub=sub, clip=clip)[
                    0
                ]
    else:
        cost[0][0] = 0

    for i in range(length + 1):
        row = cost[i]
        for k in range(states - 1):
            if row[k] < row[k + 1] or (late and row[k] == row[k + 1] < _INF):
                row[k + 1] = row[k]
                back[i][k + 1] = (i, k, False)
        for k in range(states):
            if row[k] >= _INF:
                continue
            unit = UNITS[k]
            size = len(unit)
            for step in (size, size - 1, size + 1):
                end = i + step
                if end > length:
                    continue
                total = row[k] + _unit_cost(region[i:end], unit, sub, indel)
                if total < cost[end][k]:
                    cost[end][k] = total
                    back[end][k] = (i, k, True)

    best, best_at = cost[length][states - 1], (length, states - 1)
    if open_end:
        best = _INF
        for i in range(length, length - max_edge - 1, -1):
            for k in range(states):
                tail_cost = _edge_cost(
                    region[i:], range(k, states), at_end=True, sub=sub, clip=clip
                )
                if cost[i][k] + tail_cost[0] < best:
                    best, best_at = cost[i][k] + tail_cost[0], (i, k)
    if best >= _INF:
        return None

    counts = [0] * states
    walked: list[tuple[int, int, int]] = []
    i, k = best_at
    tail = region[i:] if open_end else ""
    while (step_back := back[i][k]) is not None:
        previous_i, previous_k, consumed = step_back
        if consumed:
            counts[k] += 1
            walked.append((k, i - previous_i, cost[i][k] - cost[previous_i][previous_k]))
        i, k = previous_i, previous_k
    head = region[:i] if open_start else ""

    clipped = 0
    if head:
        cost_head, was_clip = _edge_cost(head, range(k + 1), at_end=False, sub=sub, clip=clip)
        clipped += cost_head if was_clip else 0
    if tail:
        last = best_at[1]
        cost_tail, was_clip = _edge_cost(tail, range(last, states), at_end=True, sub=sub, clip=clip)
        clipped += cost_tail if was_clip else 0
    errors = (best - clipped) / sub
    return Decomposition(tuple(counts), errors, tuple(reversed(walked)), head, tail)  # type: ignore[arg-type]


def _error_free_run(units: list[tuple[int, int, int]], needed: int) -> bool:
    run = 0
    for _, size, unit_cost in units:
        if unit_cost:
            return False
        run += size
        if run >= needed:
            return True
    return False


def field_status(
    d: Decomposition, open_start: bool, open_end: bool, confirm_bases: int = 12
) -> tuple[FieldStatus, ...]:
    """Which counts a read pins down exactly, only bounds below, or does not see.

    At an open end the last tract seen is a lower bound and later tracts are unseen.
    Two things push that frontier back towards the anchored end:

    - A zero count is only exact if what follows could not be the start of that
      tract's unit: ``...CAACAG CCG`` may be a truncated CCGCCA. In reverse at an
      open start, ``CAG CCGCCA...`` may follow the tail of a CAACAG.
    - A boundary needs ``confirm_bases`` error-free bases beyond it, otherwise the
      tract before it may simply continue.

    Fields left ambiguous by a tie in the decomposition are unobserved, so the other
    read strand decides.
    """
    n = len(UNITS)
    if not (open_start or open_end):
        return tuple(
            FieldStatus.UNOBSERVED if unsure else FieldStatus.EXACT for unsure in d.ambiguous
        )
    if not d.path:
        return (FieldStatus.UNOBSERVED,) * n

    exact: range
    bound: int | None
    if open_end:
        # Fields before `upto` are exact; `bound`, if any, is a lower bound.
        upto, bound = d.path[-1], d.path[-1]
        for z in range(d.path[-1]):
            after = "".join(UNITS[k] for k in d.path if k > z) + d.tail
            if not d.counts[z] and (after.startswith(UNITS[z]) or UNITS[z].startswith(after)):
                upto, bound = z, None
                break
        while upto > 0:
            j = upto - 1
            if _error_free_run([u for u in d.units if u[0] > j], confirm_bases):
                break
            upto, bound = j, (j if d.counts[j] else None)
        exact = range(upto)
    else:
        # Fields from `start` on are exact; `bound`, if any, is a lower bound.
        start, bound = d.path[0] + 1, d.path[0]
        for z in range(n - 1, d.path[0], -1):
            before = d.head + "".join(UNITS[k] for k in d.path if k < z)
            if not d.counts[z] and (before.endswith(UNITS[z]) or UNITS[z].endswith(before)):
                start, bound = z + 1, None
                break
        while start < n:
            j = start
            if _error_free_run([u for u in reversed(d.units) if u[0] < j], confirm_bases):
                break
            start, bound = j + 1, (j if d.counts[j] else None)
        exact = range(start, n)

    status = [FieldStatus.UNOBSERVED] * n
    for k in exact:
        status[k] = FieldStatus.EXACT
    if bound is not None:
        status[bound] = FieldStatus.LOWER_BOUND
    for k, unsure in enumerate(d.ambiguous):
        if unsure:
            status[k] = FieldStatus.UNOBSERVED
    return tuple(status)


class RepeatParser:
    """Parse forward-orientation reads into repeat-structure observations."""

    def __init__(
        self,
        amplicon: AmpliconSpec = HTT_AMPLICON,
        settings: ParserSettings | None = None,
        cache_size: int = 200_000,
    ) -> None:
        self.amplicon = amplicon
        self.settings = settings or ParserSettings()
        self._five = amplicon.five_prime_anchor(self.settings.anchor_length)
        self._three = amplicon.three_prime_anchor(self.settings.anchor_length)
        self._cache: dict[str, ReadParse] = {}
        self._cache_size = cache_size

    def parse(self, read: str) -> ReadParse:
        if (cached := self._cache.get(read)) is not None:
            return cached
        result = self._parse(read)
        if len(self._cache) >= self._cache_size:
            self._cache.clear()
        self._cache[read] = result
        return result

    def _parse(self, read: str) -> ReadParse:
        s = self.settings
        five = find_anchor(
            read,
            self._five,
            max_mismatches=s.anchor_max_mismatches,
            min_partial=s.min_partial_anchor,
            partial_side="left",
        )
        begin = five + len(self._five) if five is not None else 0
        three = find_anchor(
            read,
            self._three,
            start=begin,
            max_mismatches=s.anchor_max_mismatches,
            min_partial=s.min_partial_anchor,
            partial_side="right",
        )
        if five is None and three is None:
            return ReadParse(ReadOutcome.NO_ANCHOR)

        end = three if three is not None else len(read)
        region = read[begin:end]
        open_start, open_end = five is None, three is None
        outcome = ReadOutcome.TRUNCATED if (open_start or open_end) else ReadOutcome.SPANNING
        if not region:
            # Adjacent anchors in a spanning read cannot be an HTT allele.
            if outcome is ReadOutcome.SPANNING:
                return ReadParse(ReadOutcome.NONCONFORMING)
            return ReadParse(outcome)

        parsed = decompose(region, open_start, open_end, s)
        if parsed is None:
            return ReadParse(ReadOutcome.NONCONFORMING, region=region)
        observation = Observation(
            parsed.counts, field_status(parsed, open_start, open_end, s.confirm_bases)
        )
        if parsed.errors > s.max_errors + s.error_rate * len(region):
            outcome = ReadOutcome.NONCONFORMING
        return ReadParse(outcome, observation, parsed.errors, region)
