"""SNPs!! yeah!!

Only substitutions are called, for now; not insertions or deletions.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict

import numpy as np

from ...amplicon import AmpliconSpec
from ...counts import SampleCounts
from ...structure import FIELDS, UNITS, AlleleStructure, FieldStatus, Observation
from ...substitutions import FIVE_PRIME, SEGMENTS, THREE_PRIME, Site
from .fit import _Fit
from .likelihood import _Model, _responsibilities
from .results import Candidate, Snp
from .settings import CallerSettings
from .table import _Table

# Alleles of the same structure within this many CAG of each other can't be told apart
# by their molecules, as in the exact-N step (lengths._local_n).
_CLOSE = 3
_READ_IN_FULL = (FieldStatus.EXACT, FieldStatus.UNCONFIRMED)
_PHRED = 10 / math.log(10)


def _call_snps(
    counts: SampleCounts, table: _Table, fit: _Fit, settings: CallerSettings
) -> tuple[Snp, ...] | None:
    if not counts.has_snp_evidence:
        return None
    model = fit.model
    responsibility = _responsibilities(table, model, fit.params)
    # Each allele on its own where their molecules can be told apart, otherwise together.
    groups: list[tuple[int | None, tuple[Candidate, ...], np.ndarray]]
    if _separable(model):
        groups = [
            (index, (allele,), share)
            for index, (allele, share) in enumerate(zip(model.alleles, responsibility, strict=True))
        ]
    else:
        groups = [(None, model.alleles, np.sum(responsibility, axis=0))]

    coverage = _Coverage(counts)
    found: list[Snp] = []
    for allele, alleles, share in groups:
        weight = {o: float(share[i]) for i, o in enumerate(table.observations)}
        found.extend(_allele_snps(counts, coverage, weight, allele, alleles, settings))
    order = {segment: i for i, segment in enumerate(SEGMENTS)}
    return tuple(
        sorted(
            found,
            key=lambda s: (
                s.allele is None,
                s.allele or 0,
                order[s.site.segment],
                s.site.unit,
                s.site.offset,
            ),
        )
    )


def _separable(model: _Model) -> bool:
    if not model.heterozygous:
        return False
    first, second = model.alleles
    if first.beyond_read_length or second.beyond_read_length:
        return True
    same = first.structure.counts[1:] == second.structure.counts[1:]
    return not (same and abs(first.structure.cag - second.structure.cag) <= _CLOSE)


class _Coverage:
    """Molecules per observation, in all and with each flank reached."""

    def __init__(self, counts: SampleCounts) -> None:
        self.molecules: Counter[Observation] = Counter()
        self.five: Counter[Observation] = Counter()
        self.three: Counter[Observation] = Counter()
        for (observation, five, three), n in counts.flank_coverage.items():
            self.molecules[observation] += n
            if five:
                self.five[observation] += n
            if three:
                self.three[observation] += n

    def depth(self, segment: str, unit: int, weight: dict[Observation, float]) -> float:
        """An allele's molecules whose reads covered ``unit`` of ``segment``."""
        if segment == FIVE_PRIME:
            return sum(weight.get(o, 0.0) * n for o, n in self.five.items())
        if segment == THREE_PRIME:
            return sum(weight.get(o, 0.0) * n for o, n in self.three.items())
        tract = FIELDS.index(segment)
        return sum(
            weight.get(o, 0.0) * n
            for o, n in self.molecules.items()
            if o.status[tract] in _READ_IN_FULL and o.counts[tract] >= unit
        )


def _allele_snps(
    counts: SampleCounts,
    coverage: _Coverage,
    weight: dict[Observation, float],
    allele: int | None,
    alleles: tuple[Candidate, ...],
    settings: CallerSettings,
) -> list[Snp]:
    # Molecules showing each substitution, by (segment, offset, base) and unit.
    shown: defaultdict[tuple[str, int, str], defaultdict[int, float]] = defaultdict(
        lambda: defaultdict(float)
    )
    for (observation, site, base), n in counts.substitutions.items():
        if (w := weight.get(observation, 0.0) * n) > 0:
            shown[(site.segment, site.offset, base)][site.unit] += w

    smallest = settings.snp_partial_share * settings.snp_min_depth
    window = settings.snp_unit_window
    found = []
    for (segment, offset, base), by_unit in shown.items():
        flank = segment in (FIVE_PRIME, THREE_PRIME)
        reach = 0 if flank else window
        while by_unit:
            unit = max(by_unit, key=by_unit.__getitem__)
            molecules = sum(by_unit.pop(u, 0.0) for u in range(unit - reach, unit + reach + 1))
            if molecules < smallest:
                continue
            depth = coverage.depth(segment, max(1, unit - reach), weight)
            snp = _snp(
                Site(segment, unit, offset),
                base,
                min(molecules, depth),
                depth,
                2 * reach + 1,
                allele,
                alleles,
                settings,
            )
            if snp is not None:
                found.append(snp)
    return found


def _snp(
    site: Site,
    base: str,
    molecules: float,
    depth: float,
    units: int,
    allele: int | None,
    alleles: tuple[Candidate, ...],
    settings: CallerSettings,
) -> Snp | None:
    if depth < settings.snp_min_depth:
        return None
    share = molecules / depth
    error = min(0.5, settings.snp_error_rate * units)
    if share < settings.snp_partial_share or share <= error:
        return None
    # Likelihood ratio of "this share of the molecules carry it" against "errors", with
    # the molecules weighted down as in the genotype fit (PCR copies aren't independent).
    scale = 1.0
    if settings.effective_molecules is not None:
        scale = min(1.0, settings.effective_molecules / depth)
    n, k = depth * scale, molecules * scale
    ratio = k * math.log(share / error)
    if k < n:
        ratio += (n - k) * math.log((1 - share) / (1 - error))
    prior = math.log(settings.snp_prior) - math.log1p(-settings.snp_prior)
    quality = min(99.0, _PHRED * float(np.logaddexp(0.0, ratio + prior)))
    if quality < settings.snp_min_quality:
        return None

    places = [_place(site, a, settings.amplicon) for a in alleles]
    first = places[0]
    if first is None:
        return None
    sequence, index = first
    ref = sequence[index]
    if ref == base:
        return None
    start = index - (index - len(settings.amplicon.five_prime_flank)) % 3
    ref_codon = sequence[start : start + 3]
    alt_codon = ref_codon[: index - start] + base + ref_codon[index - start + 1 :]
    positions = {None if p is None or p[1] is None else p[1] for p in places}
    position = positions.pop() if len(positions) == 1 else None
    if position is not None:
        position += 1
    if any(a.beyond_read_length for a in alleles) and site.segment not in (FIVE_PRIME, "cag"):
        # Past a CAG tract too long to read, the place along the sequence is unknown.
        position = None
    return Snp(
        site=site,
        ref=ref,
        alt=base,
        allele=allele,
        molecules=molecules,
        depth=depth,
        quality=quality,
        partial=allele is not None and share < settings.snp_min_share,
        ref_codon=ref_codon,
        alt_codon=alt_codon,
        position=position,
    )


def _place(site: Site, allele: Candidate, amplicon: AmpliconSpec) -> tuple[str, int] | None:
    """The allele's sequence (5' flank, repeat, 3' flank) and the site's index in it, or
    None if the allele doesn't have that unit."""
    structure: AlleleStructure = allele.structure
    five, three = amplicon.five_prime_flank, amplicon.three_prime_flank
    repeat = structure.repeat_sequence()
    sequence = five + repeat + three
    if site.segment == FIVE_PRIME:
        index = len(five) - site.unit
    elif site.segment == THREE_PRIME:
        index = len(five) + len(repeat) + site.unit - 1
    else:
        tract = FIELDS.index(site.segment)
        if not 1 <= site.unit <= structure.counts[tract]:
            return None
        before = sum(len(UNITS[k]) * structure.counts[k] for k in range(tract))
        index = len(five) + before + (site.unit - 1) * len(UNITS[tract]) + site.offset
    if not 0 <= index < len(sequence):
        return None
    return sequence, index
