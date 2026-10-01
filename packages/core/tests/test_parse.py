import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from scalehd.amplicon import HTT_AMPLICON
from scalehd.parse import ParserSettings, ReadOutcome, RepeatParser, find_anchor
from scalehd.seqio import reverse_complement
from scalehd.structure import AlleleStructure, FieldStatus

FIVE = HTT_AMPLICON.five_prime_flank
THREE = HTT_AMPLICON.three_prime_flank
SPACER = "ACG"
READTHROUGH = "CTGTCTCTTATACACATCT"


def amplicon(label: str) -> str:
    return HTT_AMPLICON.sequence(AlleleStructure.from_label(label).repeat_sequence())


@pytest.fixture
def parser() -> RepeatParser:
    return RepeatParser()


@pytest.mark.parametrize(
    "label",
    [
        "17_1_1_7_2",
        "43_1_1_7_2",
        "42_0_1_7_2",
        "19_2_1_10_2",
        "20_1_0_8_2",
        "20_1_1_7_3",
        "6_1_1_7_2",
        "35_1_2_7_2",
    ],
)
def test_spanning_read_is_exact(parser: RepeatParser, label: str) -> None:
    result = parser.parse(SPACER + amplicon(label) + READTHROUGH)
    assert result.outcome is ReadOutcome.SPANNING
    assert result.observation is not None
    assert result.observation.label == label
    assert result.errors == 0


def test_substitution_in_repeat_is_tolerated(parser: RepeatParser) -> None:
    read = amplicon("30_1_1_7_2")
    i = len(FIVE) + 3 * 10 + 1  # middle base of the 11th CAG
    mutated = read[:i] + "T" + read[i + 1 :]
    result = parser.parse(mutated)
    assert result.outcome is ReadOutcome.SPANNING
    assert result.observation is not None
    assert result.observation.label == "30_1_1_7_2"
    assert result.errors == 1


def test_single_base_deletion_is_tolerated(parser: RepeatParser) -> None:
    read = amplicon("30_1_1_7_2")
    i = len(FIVE) + 3 * 12 + 1
    result = parser.parse(read[:i] + read[i + 1 :])
    assert result.observation is not None
    assert result.observation.label == "30_1_1_7_2"


def test_anchor_mismatches_tolerated(parser: RepeatParser) -> None:
    read = amplicon("20_1_1_7_2")
    i = len(FIVE) - 5  # inside the 5' anchor
    result = parser.parse(read[:i] + ("A" if read[i] != "A" else "C") + read[i + 1 :])
    assert result.outcome is ReadOutcome.SPANNING


def test_right_truncated_read_bounds_cag(parser: RepeatParser) -> None:
    read = (SPACER + amplicon("95_1_1_7_2"))[:300]
    result = parser.parse(read)
    assert result.outcome is ReadOutcome.TRUNCATED
    assert result.observation is not None
    assert result.observation.label == "83+_?_?_?_?"


def test_left_truncated_read_base_pairing_gives_exact_ccg(parser: RepeatParser) -> None:
    r2 = (SPACER + reverse_complement(amplicon("95_1_1_10_3")))[:300]
    result = parser.parse(reverse_complement(r2))
    assert result.outcome is ReadOutcome.TRUNCATED
    assert result.observation is not None
    cag, *rest = result.observation.status
    assert cag is FieldStatus.LOWER_BOUND
    assert all(s is FieldStatus.EXACT for s in rest)
    assert result.observation.counts[1:] == (1, 1, 10, 3)


def test_partial_three_prime_anchor_counts_as_spanning(parser: RepeatParser) -> None:
    read = FIVE + AlleleStructure(20).repeat_sequence() + THREE[:12]
    result = parser.parse(read)
    assert result.outcome is ReadOutcome.SPANNING
    assert result.observation is not None
    assert result.observation.label == "20_1_1_7_2"


def test_off_target_read(parser: RepeatParser) -> None:
    assert parser.parse("ACGT" * 75).outcome is ReadOutcome.NO_ANCHOR


def test_junk_between_anchors_is_nonconforming(parser: RepeatParser) -> None:
    read = FIVE + "TTTTGGGGAAAATTTTGGGGAAAATTTTGGGG" + THREE
    assert parser.parse(read).outcome is ReadOutcome.NONCONFORMING


def test_adjacent_anchors_are_nonconforming(parser: RepeatParser) -> None:
    assert parser.parse(FIVE + THREE).outcome is ReadOutcome.NONCONFORMING


@pytest.mark.parametrize(
    ("tail", "expected"),
    [
        # A trailing CCG may be the start of a CCGCCA, and the CAG boundary is only
        # nine bases from the end of the read, too few to confirm it.
        ("CAG" * 20 + "CAACAG" + "CCG", "20+_?_?_?_?"),
        # A trailing CAA is a partial CAACAG, not a miscalled CAG.
        ("CAG" * 20 + "CAA", "20+_?_?_?_?"),
        # Loss of interruption is still called; CCGCCA is only a bound because
        # just nine bases follow it.
        ("CAG" * 40 + "CCGCCA" + "CCG" * 3, "40_0_1+_?_?"),
        ("CAG" * 40 + "CCGCCA" + "CCG" * 5, "40_0_1_5+_?"),
    ],
)
def test_open_end_ambiguity(parser: RepeatParser, tail: str, expected: str) -> None:
    result = parser.parse(FIVE + tail)
    assert result.observation is not None
    assert result.observation.label == expected


def test_open_start_does_not_assume_loi(parser: RepeatParser) -> None:
    # "CAG CCGCCA ..." at a truncated start may be the tail of a CAACAG.
    result = parser.parse("CAGCCGCCA" + "CCG" * 7 + "CCT" * 2 + THREE)
    assert result.observation is not None
    assert result.observation.status[:2] == (FieldStatus.UNOBSERVED,) * 2


def test_find_anchor_partial_left() -> None:
    anchor = HTT_AMPLICON.five_prime_anchor(20)
    read = anchor[8:] + "CAGCAG"
    assert find_anchor(read, anchor, max_mismatches=2, min_partial=10, partial_side="left") == -8
    assert find_anchor(read, anchor, max_mismatches=2) is None


def test_fake_boundary_in_read_tail_is_not_trusted(parser: RepeatParser) -> None:
    # A G>A error in the last full CAG of a truncated read looks exactly like CAACAG.
    read = FIVE + "CAG" * 80 + "CAA" + "CAG" + "CA"
    result = parser.parse(read)
    assert result.observation is not None
    assert result.observation.label == "80+_?_?_?_?"


def test_tie_at_tract_boundary_is_left_to_the_read_base_pairing(parser: RepeatParser) -> None:
    # CCA in place of the last CCG is one substitution from both CCG and CCT.
    region = "CAG" * 30 + "CAACAG" + "CCGCCA" + "CCG" * 6 + "CCA" + "CCT" * 2
    result = parser.parse(FIVE + region + THREE)
    assert result.outcome is ReadOutcome.SPANNING
    assert result.observation is not None
    assert result.observation.label == "30_1_1_?_?"


def test_cache_returns_same_result(parser: RepeatParser) -> None:
    read = amplicon("20_1_1_7_2")
    assert parser.parse(read) is parser.parse(read)


structures = st.builds(
    AlleleStructure,
    cag=st.integers(1, 120),
    caacag=st.integers(0, 3),
    ccgcca=st.integers(0, 3),
    ccg=st.integers(1, 14),
    cct=st.integers(0, 4),
)


@settings(max_examples=300, deadline=None)
@given(structure=structures, spacer=st.integers(0, 6))
def test_spanning_reads_recover_structure(structure: AlleleStructure, spacer: int) -> None:
    read = "T" * spacer + HTT_AMPLICON.sequence(structure.repeat_sequence()) + READTHROUGH
    result = RepeatParser().parse(read)
    assert result.outcome is ReadOutcome.SPANNING
    assert result.observation is not None
    assert result.observation.structure() == structure


def _consistent(
    structure: AlleleStructure, result_counts: tuple[int, ...], status: tuple[FieldStatus, ...]
) -> bool:
    for truth, seen, s in zip(structure.counts, result_counts, status, strict=True):
        if s is FieldStatus.EXACT and seen != truth:
            return False
        if s is FieldStatus.LOWER_BOUND and seen > truth:
            return False
    return True


@settings(max_examples=400, deadline=None)
@given(structure=structures, data=st.data())
def test_truncated_reads_never_overclaim(structure: AlleleStructure, data: st.DataObject) -> None:
    """Without sequencing errors, truncated reads must be consistent with the truth."""
    full = HTT_AMPLICON.sequence(structure.repeat_sequence())
    repeat_start = len(FIVE)
    repeat_end = len(full) - len(THREE)
    parser = RepeatParser(settings=ParserSettings())

    cut = data.draw(st.integers(repeat_start, repeat_end), label="r1 end")
    r1 = parser.parse(full[:cut])
    if r1.observation is not None and r1.usable:
        assert _consistent(structure, r1.observation.counts, r1.observation.status), r1

    cut = data.draw(st.integers(repeat_start, repeat_end), label="r2 start")
    r2 = parser.parse(full[cut:])
    if r2.observation is not None and r2.usable:
        assert _consistent(structure, r2.observation.counts, r2.observation.status), r2
