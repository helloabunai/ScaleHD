"""Whether a call's alleles are a simulated sample's true alleles."""

import pytest
from scalehd.genotype import Candidate
from scalehd.simulate import call_matches, true_genotype
from scalehd.structure import AlleleStructure

S = AlleleStructure.from_label


def exact(label: str) -> Candidate:
    return Candidate(S(label))


def lower_bound(label: str) -> Candidate:
    return Candidate(S(label), beyond_read_length=True)


def test_true_genotype_doubles_a_homozygote_and_sorts() -> None:
    assert true_genotype([S("21_1_1_7_2")]) == (S("21_1_1_7_2"), S("21_1_1_7_2"))
    assert true_genotype([S("43_1_1_7_2"), S("17_1_1_7_2")]) == (S("17_1_1_7_2"), S("43_1_1_7_2"))


def test_true_genotype_needs_one_or_two_alleles() -> None:
    with pytest.raises(ValueError, match="one or two"):
        true_genotype([])


@pytest.mark.parametrize(
    ("called", "truth", "matches"),
    [
        ([exact("17_1_1_7_2"), exact("43_1_1_7_2")], ["17_1_1_7_2", "43_1_1_7_2"], True),
        ([exact("43_1_1_7_2"), exact("17_1_1_7_2")], ["17_1_1_7_2", "43_1_1_7_2"], True),
        ([exact("17_1_1_7_2"), exact("44_1_1_7_2")], ["17_1_1_7_2", "43_1_1_7_2"], False),
        ([exact("21_1_1_7_2"), exact("21_1_1_7_2")], ["21_1_1_7_2", "21_1_1_7_2"], True),
        ([exact("21_1_1_7_2"), exact("22_1_1_7_2")], ["21_1_1_7_2", "21_1_1_7_2"], False),
        ([exact("17_1_1_10_2"), exact("17_1_1_7_2")], ["17_1_1_7_2", "17_1_1_10_2"], True),
        ([exact("20_1_1_7_2"), lower_bound("83_1_1_7_2")], ["20_1_1_7_2", "95_1_1_7_2"], True),
        ([exact("20_1_1_7_2"), lower_bound("95_1_1_7_2")], ["20_1_1_7_2", "95_1_1_7_2"], True),
        ([exact("20_1_1_7_2"), lower_bound("96_1_1_7_2")], ["20_1_1_7_2", "95_1_1_7_2"], False),
        ([exact("20_1_1_7_2"), lower_bound("83_1_1_10_2")], ["20_1_1_7_2", "95_1_1_7_2"], False),
    ],
)
def test_call_matches(called: list[Candidate], truth: list[str], matches: bool) -> None:
    assert call_matches(called, [S(label) for label in truth]) is matches
