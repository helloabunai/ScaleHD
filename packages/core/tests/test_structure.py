import pytest
from scalehd.amplicon import HTT_AMPLICON
from scalehd.structure import AlleleStructure, FieldStatus, Observation

# Example reference from the ScaleHD 1.x docs (legacy/docs/DataAssumptions.rst).
LEGACY_17_1_1_6_2 = (
    "GCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG"
    "CAGCAGCAGCAGCAGCAGCAGCAACAGCCGCCACCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCCGCCG"
    "CAGGCACAGCCGCTGCT"
)


def test_label_round_trip() -> None:
    allele = AlleleStructure.from_label("42_0_1_7_2")
    assert allele == AlleleStructure(42, 0, 1, 7, 2)
    assert allele.label == "42_0_1_7_2"
    assert str(allele) == "42_0_1_7_2"


@pytest.mark.parametrize("label", ["17_1_1_7", "17_1_1_7_2_0", "17_a_1_7_2", "", "-1_1_1_7_2"])
def test_bad_labels(label: str) -> None:
    with pytest.raises(ValueError, match="label"):
        AlleleStructure.from_label(label)


def test_negative_counts_rejected() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        AlleleStructure(-1)


def test_defaults_are_the_common_allele() -> None:
    allele = AlleleStructure(20)
    assert allele.label == "20_1_1_7_2"
    assert allele.is_typical


@pytest.mark.parametrize(
    ("label", "typical"),
    [
        ("20_1_1_7_2", True),
        ("20_1_1_10_2", True),
        ("42_0_1_7_2", False),
        ("19_2_1_7_2", False),
        ("20_1_0_8_2", False),
        ("20_1_1_7_3", False),
    ],
)
def test_is_typical(label: str, typical: bool) -> None:
    assert AlleleStructure.from_label(label).is_typical is typical


def test_polyglutamine_length_counts_caa() -> None:
    # Loss of the CAA interruption lengthens the pure CAG tract without changing polyQ.
    assert AlleleStructure.from_label("40_1_1_7_2").polyglutamine_length == 42
    assert AlleleStructure.from_label("42_0_1_7_2").polyglutamine_length == 42
    assert AlleleStructure.from_label("40_2_1_7_2").polyglutamine_length == 44


def test_reproduces_legacy_reference_sequence() -> None:
    allele = AlleleStructure.from_label("17_1_1_6_2")
    assert HTT_AMPLICON.sequence(allele.repeat_sequence()) == LEGACY_17_1_1_6_2


def test_with_counts() -> None:
    assert AlleleStructure(20).with_counts(cag=21, ccg=10).label == "21_1_1_10_2"


def test_observation_labels() -> None:
    exact = Observation.exact(AlleleStructure(20))
    assert exact.is_complete
    assert exact.label == "20_1_1_7_2"
    assert exact.structure() == AlleleStructure(20)

    partial = Observation(
        (84, 0, 0, 0, 0),
        (FieldStatus.LOWER_BOUND,) + (FieldStatus.UNOBSERVED,) * 4,
    )
    assert not partial.is_complete
    assert partial.label == "84+_?_?_?_?"
    with pytest.raises(ValueError, match="incomplete"):
        partial.structure()
