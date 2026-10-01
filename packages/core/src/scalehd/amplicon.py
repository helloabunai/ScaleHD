"""Amplicon definitions: the sequence either side of the repeat region."""

from __future__ import annotations

from dataclasses import dataclass

_BASES = frozenset("ACGT")


@dataclass(frozen=True, slots=True)
class AmpliconSpec:
    """Flanks of the repeat region, in forward (CAG-strand) orientation.

    Reads are located by short anchors taken from the flank ends nearest the
    repeat, so any primer pair whose product contains those anchors works.
    """

    name: str
    five_prime_flank: str
    three_prime_flank: str

    def __post_init__(self) -> None:
        for flank in (self.five_prime_flank, self.three_prime_flank):
            if not flank or not set(flank) <= _BASES:
                raise ValueError("flanks must be non-empty and contain only A, C, G, T")

    def five_prime_anchor(self, length: int) -> str:
        return self.five_prime_flank[-length:]

    def three_prime_anchor(self, length: int) -> str:
        return self.three_prime_flank[:length]

    def sequence(self, repeat: str) -> str:
        return self.five_prime_flank + repeat + self.three_prime_flank


# Flanking sequenceof every sequence in the ScaleHD 1.x reference library (legacy/ScaleHD/config).
# The GeM-HD MiSeq primers (Cell 2019; ATGAAGGCCTTCGAGTCCC / GGCTGAGGAAGCTGAGGA)
HTT_AMPLICON = AmpliconSpec(
    name="HTT exon 1",
    five_prime_flank="GCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTC",
    three_prime_flank="CAGCTTCCTCAGCCGCCGCCGCAGGCACAGCCGCTGCT",
)
