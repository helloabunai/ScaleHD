"""Amplicon definitions. The sequence either side of the repeat region, and the primers."""

from __future__ import annotations

from dataclasses import dataclass, field, replace

from .seqio import reverse_complement

_BASES = frozenset("ACGT")
# Shortest stretch of a primer that has to match the flank to place it.
_MIN_PRIMER_MATCH = 10


@dataclass(frozen=True, slots=True)
class PrimerSites:
    """Where the primers sit, as flank positions (0-based, half open).

    A primer can run into the repeat region: its last bases are then repeat bases, counted
    by ``into_repeat``.
    """

    # 5' flank positions the forward primer covers, and how far it runs into the repeat.
    forward: tuple[int, int]
    forward_into_repeat: int
    # 3' flank positions the reverse primer covers (forward strand), and how far it runs
    # back into the repeat.
    reverse: tuple[int, int]
    reverse_into_repeat: int


@dataclass(frozen=True, slots=True)
class AmpliconSpec:
    """Flanks of the repeat region, in forward (CAG-strand) orientation.

    Reads are located by short anchors taken from the flank ends nearest the
    repeat, so any primer pair whose product contains those anchors works.

    The primers are written 5' to 3' as ordered: the forward primer as on the CAG strand,
    the reverse primer as on the other strand. Bases a read gets from a primer show the
    primer's sequence rather than the sample's, so SNPs aren't called in them, and the
    anchors are kept inside the product. Without primers, the whole flanks are taken
    as the product.
    """

    name: str
    five_prime_flank: str
    three_prime_flank: str
    forward_primer: str | None = None
    reverse_primer: str | None = None
    primer_sites: PrimerSites | None = field(init=False, default=None, compare=False)

    def __post_init__(self) -> None:
        for flank in (self.five_prime_flank, self.three_prime_flank):
            if not flank or not set(flank) <= _BASES:
                raise ValueError("flanks must be non-empty and contain only A, C, G, T")
        if (self.forward_primer is None) != (self.reverse_primer is None):
            raise ValueError("give both primers, or neither")
        if self.forward_primer is not None and self.reverse_primer is not None:
            sites = _place_primers(self, self.forward_primer, self.reverse_primer)
            object.__setattr__(self, "primer_sites", sites)

    def with_primers(self, forward: str | None, reverse: str | None) -> AmpliconSpec:
        return replace(self, forward_primer=forward, reverse_primer=reverse)

    @property
    def product_flanks(self) -> tuple[str, str]:
        """The parts of the flanks inside the PCR product: from the forward primer's start
        to the reverse primer's end."""
        if self.primer_sites is None:
            return self.five_prime_flank, self.three_prime_flank
        start = self.primer_sites.forward[0]
        end = self.primer_sites.reverse[1]
        return self.five_prime_flank[start:], self.three_prime_flank[:end]

    def five_prime_anchor(self, length: int) -> str:
        return self.product_flanks[0][-length:]

    def three_prime_anchor(self, length: int) -> str:
        return self.product_flanks[1][:length]

    def sequence(self, repeat: str) -> str:
        """The PCR product for an allele's repeat sequence."""
        five, three = self.product_flanks
        return five + repeat + three

    def five_prime_window(self, trim_primers: bool, anchor_length: int) -> tuple[int, int]:
        """Where SNPs can be called on the 5' side."""
        if self.primer_sites is None:
            return max(0, len(self.five_prime_flank) - anchor_length), 0
        start, end = self.primer_sites.forward
        if not trim_primers:
            return start, 0
        return min(end, len(self.five_prime_flank)), self.primer_sites.forward_into_repeat

    def three_prime_window(self, trim_primers: bool, anchor_length: int) -> tuple[int, int]:
        """Where SNPs can be called on the 3' side."""
        if self.primer_sites is None:
            return min(anchor_length, len(self.three_prime_flank)), 0
        start, end = self.primer_sites.reverse
        if not trim_primers:
            return end, 0
        return max(start, 0), self.primer_sites.reverse_into_repeat


def _place_primers(amplicon: AmpliconSpec, forward: str, reverse: str) -> PrimerSites:
    five, three = amplicon.five_prime_flank, amplicon.three_prime_flank
    for name, primer in (("forward", forward), ("reverse", reverse)):
        if len(primer) < _MIN_PRIMER_MATCH or not set(primer) <= _BASES:
            raise ValueError(
                f"the {name} primer must be at least {_MIN_PRIMER_MATCH} bases of A, C, G, T"
            )

    # The forward primer sits in the 5' flank, or runs from it into the repeat region.
    forward_site: tuple[int, int] | None = None
    forward_into = 0
    if (at := five.find(forward)) >= 0:
        forward_site = (at, at + len(forward))
    else:
        for into in range(1, len(forward) - _MIN_PRIMER_MATCH + 1):
            if five.endswith(forward[: len(forward) - into]):
                forward_site, forward_into = (len(five) - len(forward) + into, len(five)), into
                break
    if forward_site is None:
        raise ValueError(f"the forward primer isn't in the {amplicon.name} 5' flank")

    # The reverse primer binds the other strand i.e. its reverse complement is in the 3' flank,
    # or runs back into the repeat region.
    site = reverse_complement(reverse)
    reverse_site: tuple[int, int] | None = None
    reverse_into = 0
    if (at := three.find(site)) >= 0:
        reverse_site = (at, at + len(site))
    else:
        for into in range(1, len(site) - _MIN_PRIMER_MATCH + 1):
            if three.startswith(site[into:]):
                reverse_site, reverse_into = (0, len(site) - into), into
                break
    if reverse_site is None:
        raise ValueError(f"the reverse primer isn't in the {amplicon.name} 3' flank")
    return PrimerSites(forward_site, forward_into, reverse_site, reverse_into)


# Flanking sequence of every sequence in the ScaleHD 1.x reference library (legacy/ScaleHD/config).
HTT_AMPLICON = AmpliconSpec(
    name="HTT exon 1",
    five_prime_flank="GCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTC",
    three_prime_flank="CAGCTTCCTCAGCCGCCGCCGCAGGCACAGCCGCTGCT",
)

# The GeM-HD MiSeq primers (Cell 2019). The forward primer ends 11 bases before the CAG
# tract, and the reverse primer starts 4 bases into the CCT tract, so the product holds
# no 3' flank beyond the primer.
GEM_HD_PRIMERS = ("ATGAAGGCCTTCGAGTCCC", "GGCTGAGGAAGCTGAGGA")
HTT_GEM_HD = HTT_AMPLICON.with_primers(*GEM_HD_PRIMERS)
