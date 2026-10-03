"""Which alleles to try: common structures, their neighbours, and one beyond read length."""

from __future__ import annotations

from collections import Counter

from ...counts import SampleCounts
from ...structure import AlleleStructure, FieldStatus
from .results import Candidate
from .settings import CallerSettings


def candidate_alleles(
    counts: SampleCounts, settings: CallerSettings | None = None
) -> list[Candidate]:
    settings = settings or CallerSettings()
    found: dict[Candidate, None] = {}
    common = [s for s, _ in _counted_structures(counts).most_common(settings.max_candidates)]
    for structure in common:
        found[Candidate(structure)] = None
    for structure in common[: settings.neighbour_candidates]:
        for step in (-1, 1):
            if structure.cag + step >= 1:
                found[Candidate(structure.with_counts(cag=structure.cag + step))] = None
    if (long_allele := _long_candidate(counts, settings)) is not None:
        # Lengths from where reads stop up are the beyond-read-length candidate's alone.
        # W/ an exact call there would be a single point of its range, unfairly weighted.
        limit = long_allele.structure.cag
        found = {c: None for c in found if c.structure.cag < limit}
        found[long_allele] = None
    return list(found)


def _counted_structures(counts: SampleCounts) -> Counter[AlleleStructure]:
    """Complete molecules, and those whose CAG end was seen but unconfirmed.

    Near read length most molecules at an allele's own length are the second kind.
    """
    found = Counter(counts.complete)
    for observation, n in counts.partial.items():
        cag, *rest = observation.status
        if cag is FieldStatus.UNCONFIRMED and all(x is FieldStatus.EXACT for x in rest):
            found[AlleleStructure.from_counts(observation.counts)] += n
    return found


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
