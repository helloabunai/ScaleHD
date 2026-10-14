"""SNP calling without alignment: primers, substitutions in reads, and calls per allele."""

from functools import cache
from pathlib import Path

import pytest
from scalehd.amplicon import GEM_HD_PRIMERS, HTT_AMPLICON, HTT_GEM_HD
from scalehd.counts import SampleCounts, count_reads
from scalehd.genotype import Flag, GenotypeCall, call_genotype
from scalehd.parse import ParserSettings, RepeatParser
from scalehd.seqio import FastqRecord, reverse_complement
from scalehd.simulate import SimAllele, SimulationSpec, simulate
from scalehd.structure import AlleleStructure, FieldStatus, Observation
from scalehd.substitutions import FIVE_PRIME, ReadEvidence, Site, amino_acid, molecule_substitutions

READTHROUGH = "CTGTCTCTTATACACATCT"


def _read(structure: AlleleStructure, changes: dict[int, str] | None = None) -> str:
    """A spanning GeM-HD read of ``structure``, with bases changed at product positions."""
    product = list(HTT_GEM_HD.sequence(structure.repeat_sequence()))
    for index, base in (changes or {}).items():
        product[index] = base
    return "ACG" + "".join(product) + READTHROUGH


def _cag_base(unit: int, offset: int) -> int:
    return len(HTT_GEM_HD.product_flanks[0]) + 3 * (unit - 1) + offset


def _subs(read: str, *, trim: bool = True) -> list[tuple[str, str]]:
    parser = RepeatParser(HTT_GEM_HD, ParserSettings(trim_primers=trim))
    return [(s.site.label + f".{s.site.offset}", s.base) for s in parser.parse(read).substitutions]


def test_the_gem_hd_primers_leave_11_bases_before_the_cag_and_none_after_the_repeats() -> None:
    sites = HTT_GEM_HD.primer_sites
    assert sites is not None
    assert (sites.forward, sites.forward_into_repeat) == ((18, 37), 0)
    assert (sites.reverse, sites.reverse_into_repeat) == ((0, 14), 4)
    assert HTT_GEM_HD.five_prime_window(True, 20) == (37, 0)
    assert HTT_GEM_HD.three_prime_window(True, 20) == (0, 4)
    # Untrimmed, the primer bases count too.
    assert HTT_GEM_HD.five_prime_window(False, 20) == (18, 0)
    assert HTT_GEM_HD.three_prime_window(False, 20) == (14, 0)


def test_anchors_stay_inside_the_pcr_product() -> None:
    assert HTT_GEM_HD.three_prime_anchor(20) == "CAGCTTCCTCAGCC"
    assert HTT_GEM_HD.five_prime_anchor(20) == HTT_AMPLICON.five_prime_anchor(20)
    product = HTT_GEM_HD.sequence(AlleleStructure(20).repeat_sequence())
    assert product.startswith(GEM_HD_PRIMERS[0])
    assert product.endswith(reverse_complement(GEM_HD_PRIMERS[1]))


@pytest.mark.parametrize(
    ("forward", "reverse", "message"),
    [
        ("ACGTACGTACGTACGT", GEM_HD_PRIMERS[1], "forward primer isn't"),
        (GEM_HD_PRIMERS[0], "ACGTACGTACGTACGT", "reverse primer isn't"),
        ("ATGAAG", GEM_HD_PRIMERS[1], "at least 10 bases"),
    ],
)
def test_primers_must_be_found_in_the_flanks(forward: str, reverse: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        HTT_AMPLICON.with_primers(forward, reverse)


def test_a_read_shows_substitutions_where_its_structure_says_otherwise() -> None:
    structure = AlleleStructure(20)
    assert _subs(_read(structure)) == []
    assert _subs(_read(structure, {_cag_base(12, 2): "T"})) == [("CAG unit 12.2", "T")]
    # The 5' flank's last bases are outside the primer; the primer's are left out.
    five = len(HTT_GEM_HD.product_flanks[0])
    assert _subs(_read(structure, {five - 3: "A"})) == [("5' flank -3.0", "A")]
    assert _subs(_read(structure, {5: "A"})) == []
    assert _subs(_read(structure, {5: "A"}), trim=False) == [("5' flank -25.0", "A")]
    # The reverse primer covers the last 4 bases of the repeat region.
    last = five + len(structure.repeat_sequence()) - 2
    assert _subs(_read(structure, {last: "A"})) == []


def test_a_read_that_starts_inside_the_cag_tract_numbers_it_from_the_end() -> None:
    structure = AlleleStructure(90)
    read = _read(structure, {_cag_base(88, 2): "T"})[3:]
    tail = reverse_complement(reverse_complement(read)[:150])
    assert _subs(tail) == [("CAG unit -3.2", "T")]


def test_a_molecule_keeps_a_substitution_its_reads_agree_on_or_one_good_read_shows() -> None:
    structure = AlleleStructure(20)
    parser = RepeatParser(HTT_GEM_HD)
    snp, plain = _read(structure, {_cag_base(5, 1): "T"}), _read(structure)
    observation = Observation.exact(structure)
    good, fair, poor = "I" * len(snp), "5" * len(snp), "+" * len(snp)

    def kept(*reads: tuple[str, str]) -> frozenset[tuple[Site, str]]:
        evidence = [ReadEvidence(read, quality, parser.parse(read)) for read, quality in reads]
        return molecule_substitutions(evidence, observation)

    shown = frozenset({(Site("cag", 5, 1), "T")})
    # Two fair (Q20) reads that agree, or one good (Q40) read alone.
    assert kept((snp, fair), (snp, fair)) == shown
    assert kept((snp, good)) == shown
    assert kept((snp, fair)) == frozenset()
    assert kept((snp, poor), (snp, poor)) == frozenset()
    # A read showing the reference base outvotes it.
    assert kept((snp, good), (plain, good)) == frozenset()


def test_snp_evidence_round_trips_through_counts_files(tmp_path: Path) -> None:
    counts = SampleCounts()
    exact = Observation.exact(AlleleStructure(20))
    partial = Observation((60, 1, 1, 7, 2), (FieldStatus.LOWER_BOUND,) + (FieldStatus.EXACT,) * 4)
    counts.add(exact, substitutions=[(Site("cag", 5, 1), "T")], flanks=(True, True))
    counts.add(partial, substitutions=[(Site(FIVE_PRIME, 3), "A")], flanks=(True, False))
    path = tmp_path / "s.counts.json"
    counts.write_json(path)
    again = SampleCounts.read_json(path)
    assert again == counts
    assert again.has_snp_evidence


@cache
def _called(*alleles: SimAllele, pairs: int = 1500) -> GenotypeCall:
    spec = SimulationSpec(alleles, pairs=pairs, seed=4, amplicon=HTT_GEM_HD, pcr_error=1e-4)
    sample = simulate(spec)
    pairs_read = zip(sample.r1, sample.r2, strict=True)
    return call_genotype(count_reads(pairs_read, RepeatParser(HTT_GEM_HD)))


def test_each_allele_gets_its_own_snps() -> None:
    call = _called(
        SimAllele(AlleleStructure(17), snps=((Site(FIVE_PRIME, 3), "A"),)),
        SimAllele(AlleleStructure(43), snps=((Site("cag", 12, 2), "T"),)),
    )
    assert call.snps is not None
    found = [(s.allele, s.label, s.legacy_label, s.partial) for s in call.snps]
    assert found == [
        (0, "5' flank -3: TTC>ATC, Phe>Ile", "T->A@46", False),
        (1, "CAG unit 12: CAG>CAT, Gln>His", "G->T@84", False),
    ]
    assert all(s.quality >= 20 and s.share > 0.9 for s in call.snps)
    assert Flag.SNP in call.flags
    assert Flag.PARTIAL_SNP not in call.flags
    assert call.to_dict()["snps"][1]["alt_amino_acid"] == "His"


def test_no_snps_are_called_from_copying_and_reading_errors() -> None:
    call = _called(SimAllele(AlleleStructure(17)), SimAllele(AlleleStructure(43)))
    assert call.snps == ()
    assert Flag.SNP not in call.flags


def test_a_snp_on_one_of_two_alleles_that_cant_be_told_apart_has_no_allele() -> None:
    call = _called(
        SimAllele(AlleleStructure(20)),
        SimAllele(AlleleStructure(20), snps=((Site("ccg", 3, 0), "T"),)),
    )
    assert call.snps is not None
    (snp,) = call.snps
    assert (snp.allele, snp.label, snp.partial) == (None, "CCG unit 3: CCG>TCG, Pro>Ser", False)
    assert 0.3 < snp.share < 0.7


def test_counts_without_snp_evidence_give_no_snp_call() -> None:
    counts = SampleCounts()
    for _ in range(600):
        counts.add(Observation.exact(AlleleStructure(17)))
        counts.add(Observation.exact(AlleleStructure(21)))
    assert call_genotype(counts).snps is None
    assert call_genotype(counts).to_dict()["snps"] is None


def test_codons_translate_with_the_standard_code() -> None:
    assert [amino_acid(c) for c in ("CAG", "CAA", "CCG", "ATG", "TAA", "TGG")] == [
        "Gln",
        "Gln",
        "Pro",
        "Met",
        "Ter",
        "Trp",
    ]


def test_records_carry_base_qualities_into_counting() -> None:
    structure = AlleleStructure(20)
    snp = _read(structure, {_cag_base(5, 1): "T"})
    pair = (FastqRecord("r", snp, "+" * len(snp)), None)
    counts = count_reads([pair], RepeatParser(HTT_GEM_HD))
    assert not counts.substitutions
    pair = (FastqRecord("r", snp, "I" * len(snp)), None)
    counts = count_reads([pair], RepeatParser(HTT_GEM_HD))
    assert list(counts.substitutions) == [(Observation.exact(structure), Site("cag", 5, 1), "T")]


def test_a_snp_only_some_of_an_alleles_molecules_show_is_partial() -> None:
    snp = ((Site("cag", 12, 2), "T"),)
    call = _called(
        SimAllele(AlleleStructure(17)),
        SimAllele(AlleleStructure(43), 0.7),
        SimAllele(AlleleStructure(43), 0.3, snps=snp),
    )
    assert call.snps is not None
    (found,) = call.snps
    assert (found.allele, found.partial) == (1, True)
    assert 0.15 < found.share < 0.45
    assert Flag.PARTIAL_SNP in call.flags
