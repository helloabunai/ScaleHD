"""Single-base substitutions in reads, placed by repeat structure.

There is no alignment. The parser has already placed each read by its anchors and split
the repeat region unit by unit, so every base's place in the allele is known. A
substitution is a base that differs from what the read's own structure says is there:

- in a flank, a base that differs from the amplicon's flank,
- in the repeat region, a base that differs from its unit, in a unit read at its full
  length (CAG read as CAT).

Only substitutions are supported for now
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

from .structure import FIELDS, UNITS, Counts, FieldStatus, Observation

if TYPE_CHECKING:
    from .amplicon import AmpliconSpec
    from .parse import Decomposition, ReadParse

FIVE_PRIME = "five_prime"
THREE_PRIME = "three_prime"
SEGMENTS: tuple[str, ...] = (FIVE_PRIME, *FIELDS, THREE_PRIME)
SEGMENT_NAMES = {FIVE_PRIME: "5' flank", THREE_PRIME: "3' flank"} | dict(
    zip(FIELDS, UNITS, strict=True)
)

# A read whose flank differs at more places than this is misplaced (an insertion or
# deletion shifts every base after it), so its flank substitutions are left out.
_MAX_FLANK_SUBSTITUTIONS = 3
# Statuses of a tract that some read saw from end to end.
_READ_IN_FULL = (FieldStatus.EXACT, FieldStatus.UNCONFIRMED)
# Base quality assumed for reads given without qualities.
_UNKNOWN_QUALITY = 40


@dataclass(frozen=True, slots=True, order=True)
class Site:
    """A base in an allele, placed by the repeat structure rather than by its number
    along the sequence, so it keeps its place whatever the tract lengths.
    """

    segment: str
    unit: int
    offset: int = 0

    @property
    def in_flank(self) -> bool:
        return self.segment in (FIVE_PRIME, THREE_PRIME)

    @property
    def label(self) -> str:
        if self.segment == FIVE_PRIME:
            return f"5' flank -{self.unit}"
        if self.segment == THREE_PRIME:
            return f"3' flank +{self.unit}"
        return f"{SEGMENT_NAMES[self.segment]} unit {self.unit}"

    @property
    def key(self) -> str:
        return f"{self.segment}:{self.unit}:{self.offset}"

    @classmethod
    def from_key(cls, key: str) -> Site:
        segment, unit, offset = key.split(":")
        if segment not in SEGMENTS:
            raise ValueError(f"unknown segment in site {key!r}")
        return cls(segment, int(unit), int(offset))


@dataclass(frozen=True, slots=True)
class Substitution:
    site: Site
    base: str
    # Where it is in the read (forward orientation), for the base's quality.
    position: int


@dataclass(frozen=True, slots=True)
class Windows:
    """Where SNPs can be called in reads, from AmpliconSpec.five_prime_window and
    three_prime_window.
    """

    five_prime_flank: str
    three_prime_flank: str
    five_start: int
    five_masked: int
    three_end: int
    three_masked: int

    @classmethod
    def of(cls, amplicon: AmpliconSpec, *, trim_primers: bool, anchor_length: int) -> Windows:
        five_start, five_masked = amplicon.five_prime_window(trim_primers, anchor_length)
        three_end, three_masked = amplicon.three_prime_window(trim_primers, anchor_length)
        return cls(
            amplicon.five_prime_flank,
            amplicon.three_prime_flank,
            five_start,
            five_masked,
            three_end,
            three_masked,
        )


def read_substitutions(
    read: str,
    begin: int,
    end: int,
    five_anchored: bool,
    three_anchored: bool,
    decomposition: Decomposition | None,
    windows: Windows,
) -> tuple[Substitution, ...]:
    """Every substitution in a parsed read, where SNPs can be called. Outside the primers
    (unless they aren't trimmed) and inside the PCR product."""
    found: list[Substitution] = []
    if five_anchored:
        flank = windows.five_prime_flank
        first = begin - len(flank)
        low = max(windows.five_start, -first)
        if read[first + low : begin] != flank[low:]:
            here = [
                Substitution(Site(FIVE_PRIME, len(flank) - f), read[first + f], first + f)
                for f in range(low, len(flank))
                if read[first + f] != flank[f] and read[first + f] != "N"
            ]
            if len(here) <= _MAX_FLANK_SUBSTITUTIONS:
                found.extend(here)
    if three_anchored and windows.three_end:
        flank = windows.three_prime_flank
        high = min(windows.three_end, len(read) - end)
        if read[end : end + high] != flank[:high]:
            here = [
                Substitution(Site(THREE_PRIME, f + 1), read[end + f], end + f)
                for f in range(high)
                if read[end + f] != flank[f] and read[end + f] != "N"
            ]
            if len(here) <= _MAX_FLANK_SUBSTITUTIONS:
                found.extend(here)

    if decomposition is not None and decomposition.errors > 0:
        # Repeat bases that come from a primer running into the repeat region.
        lowest = begin + windows.five_masked if five_anchored else begin
        highest = end - windows.three_masked if three_anchored else end
        found.extend(
            _repeat_substitutions(read, begin, five_anchored, decomposition, lowest, highest)
        )
    return tuple(found)


def _repeat_substitutions(
    read: str, begin: int, five_anchored: bool, d: Decomposition, lowest: int, highest: int
) -> list[Substitution]:
    seen = [0] * len(UNITS)
    raw: list[tuple[int, int, int, str, int]] = []
    position = begin + len(d.head)
    for tract, size, cost in d.units:
        index = seen[tract]
        seen[tract] += 1
        unit = UNITS[tract]
        if cost and size == len(unit) and not d.ambiguous[tract]:
            for offset, (base, expected) in enumerate(
                zip(read[position : position + size], unit, strict=True)
            ):
                at = position + offset
                if base != expected and base != "N" and lowest <= at < highest:
                    raw.append((tract, index, offset, base, at))
        position += size
    # A read that starts inside the repeat region didn't see its first tract's start, so
    # that tract is numbered from its end.
    first = d.units[0][0] if d.units else None
    found = []
    for tract, index, offset, base, at in raw:
        from_end = not five_anchored and tract == first
        unit_number = index - seen[tract] if from_end else index + 1
        found.append(Substitution(Site(FIELDS[tract], unit_number, offset), base, at))
    return found


def site_position(parse: ReadParse, length: int, site: Site, tract_counts: Counts) -> int | None:
    """Where ``site`` is in a parsed read of ``length`` bases, or None if the read doesn't
    show it. ``site`` is numbered from its tract's start so ``tract_counts`` gives the tract
    lengths, for a read that didn't see a tract's start."""
    if site.segment == FIVE_PRIME:
        at = parse.begin - site.unit
        return at if parse.five_anchored and 0 <= at < length else None
    if site.segment == THREE_PRIME:
        at = parse.end + site.unit - 1
        return at if parse.three_anchored and 0 <= at < length else None
    d = parse.decomposition
    if d is None:
        return None
    tract = FIELDS.index(site.segment)
    if d.ambiguous[tract]:
        return None
    # Units come in tract order, so this tract's are a run of d.counts[tract] of them.
    in_read = d.counts[tract]
    from_end = not parse.five_anchored and bool(d.units) and d.units[0][0] == tract
    index = in_read - 1 - (tract_counts[tract] - site.unit) if from_end else site.unit - 1
    if not 0 <= index < in_read:
        return None
    unit = sum(d.counts[:tract]) + index
    size = d.units[unit][1]
    if size != len(UNITS[tract]):
        return None
    at = parse.begin + len(d.head) + sum(step for _, step, _ in d.units[:unit])
    return at + site.offset


@dataclass(frozen=True, slots=True)
class ReadEvidence:
    """One read of a molecule, parsed, with its base qualities (Phred+33) if known."""

    read: str
    quality: str | None
    parse: ReadParse


def molecule_substitutions(
    reads: Sequence[ReadEvidence],
    observation: Observation,
    *,
    min_quality: int = 30,
    min_pair_quality: int = 20,
) -> frozenset[tuple[Site, str]]:
    """The substitutions a molecule's reads stand by, numbered from each tract's start."""
    counts = observation.counts

    def placed(site: Site) -> Site | None:
        if site.in_flank:
            return site
        tract = FIELDS.index(site.segment)
        if observation.status[tract] not in _READ_IN_FULL:
            return None
        unit = site.unit if site.unit > 0 else counts[tract] + site.unit + 1
        return Site(site.segment, unit, site.offset) if 1 <= unit <= counts[tract] else None

    def quality(evidence: ReadEvidence, at: int) -> int:
        return _UNKNOWN_QUALITY if evidence.quality is None else ord(evidence.quality[at]) - 33

    kept: set[tuple[Site, str]] = set()
    for i, evidence in enumerate(reads):
        for substitution in evidence.parse.substitutions:
            site = placed(substitution.site)
            own = quality(evidence, substitution.position)
            if site is None or own < min_pair_quality:
                continue
            confirmed = contradicted = False
            for j, other in enumerate(reads):
                if j == i:
                    continue
                at = site_position(other.parse, len(other.read), site, counts)
                if at is None or quality(other, at) < min_pair_quality:
                    continue
                if other.read[at] == substitution.base:
                    confirmed = True
                else:
                    contradicted = True
            if not contradicted and (confirmed or own >= min_quality):
                kept.add((site, substitution.base))
    return frozenset(kept)


_CODON_BASES = "TCAG"
_ONE_LETTER = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_THREE_LETTER = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
    "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
    "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter",
}  # fmt: skip
_CODONS = {
    a + b + c: _THREE_LETTER[_ONE_LETTER[16 * i + 4 * j + k]]
    for i, a in enumerate(_CODON_BASES)
    for j, b in enumerate(_CODON_BASES)
    for k, c in enumerate(_CODON_BASES)
}


def amino_acid(codon: str) -> str:
    """The standard genetic code's amino acid for a codon, as three letters (Ter for stop)."""
    return _CODONS[codon]
