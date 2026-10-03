from scalehd.pairs import DiscordancePolicy, join_read_base_pairings
from scalehd.structure import AlleleStructure, FieldStatus, Observation

E, L, U = FieldStatus.EXACT, FieldStatus.LOWER_BOUND, FieldStatus.UNOBSERVED
UC = FieldStatus.UNCONFIRMED


def exact(label: str) -> Observation:
    return Observation.exact(AlleleStructure.from_label(label))


def test_agreeing_read_base_pairings() -> None:
    result = join_read_base_pairings(exact("43_1_1_7_2"), exact("43_1_1_7_2"))
    assert result.observation == exact("43_1_1_7_2")
    assert not result.is_discordant


def test_disagreeing_read_base_pairings_dropped_by_default() -> None:
    result = join_read_base_pairings(exact("43_1_1_7_2"), exact("41_2_1_7_2"))
    assert result.observation is None
    assert result.discordant == (True, True, False, False, False)


def test_prefer_policy_takes_preferred_read_base_pairing_per_field() -> None:
    result = join_read_base_pairings(
        exact("43_1_1_7_2"), exact("43_1_1_8_2"), policy=DiscordancePolicy.PREFER
    )
    assert result.observation is not None
    assert result.observation.label == "43_1_1_8_2"  # CCG comes from R2


def test_long_allele_completed_from_both_read_base_pairings() -> None:
    r1 = Observation((75, 1, 1, 2, 0), (E, E, E, L, U))
    r2 = Observation((60, 1, 1, 7, 2), (L, E, E, E, E))
    result = join_read_base_pairings(r1, r2)
    assert result.observation is not None
    assert result.observation.label == "75_1_1_7_2"
    assert not result.is_discordant


def test_lower_bound_when_neither_read_base_pairing_spans_cag() -> None:
    r1 = Observation((84, 0, 0, 0, 0), (L, U, U, U, U))
    r2 = Observation((70, 1, 1, 7, 2), (L, E, E, E, E))
    result = join_read_base_pairings(r1, r2)
    assert result.observation is not None
    assert result.observation.label == "84+_1_1_7_2"


def test_lower_bound_above_exact_is_discordant() -> None:
    r1 = Observation((40, 1, 1, 7, 2), (E, E, E, E, E))
    r2 = Observation((45, 1, 1, 7, 2), (L, E, E, E, E))
    assert join_read_base_pairings(r1, r2).observation is None


def test_single_read_base_pairing() -> None:
    assert join_read_base_pairings(exact("20_1_1_7_2"), None).observation == exact("20_1_1_7_2")
    assert join_read_base_pairings(None, None).observation is None


def test_unconfirmed_count_kept_when_the_other_read_does_not_see_its_end() -> None:
    r1 = Observation((80, 0, 0, 0, 0), (UC, U, U, U, U))
    r2 = Observation((70, 1, 1, 7, 2), (L, E, E, E, E))
    result = join_read_base_pairings(r1, r2)
    assert result.observation is not None
    assert result.observation.label == "80~_1_1_7_2"
    assert not result.is_discordant


def test_unconfirmed_count_agreeing_with_an_exact_one_is_exact() -> None:
    r1 = Observation((80, 0, 0, 0, 0), (UC, U, U, U, U))
    result = join_read_base_pairings(r1, exact("80_1_1_7_2"))
    assert result.observation == exact("80_1_1_7_2")


def test_unconfirmed_count_below_an_exact_one_was_a_sequencing_error() -> None:
    r1 = Observation((78, 0, 0, 0, 0), (UC, U, U, U, U))
    result = join_read_base_pairings(r1, exact("80_1_1_7_2"))
    assert result.observation == exact("80_1_1_7_2")
    assert not result.is_discordant


def test_unconfirmed_count_above_an_exact_one_is_discordant() -> None:
    r1 = Observation((82, 0, 0, 0, 0), (UC, U, U, U, U))
    assert join_read_base_pairings(r1, exact("80_1_1_7_2")).observation is None


def test_longer_lower_bound_overrules_an_unconfirmed_count() -> None:
    # The other read saw more CAG than the unconfirmed end allows, so that end was an error.
    r1 = Observation((80, 0, 0, 0, 0), (UC, U, U, U, U))
    r2 = Observation((84, 1, 1, 7, 2), (L, E, E, E, E))
    result = join_read_base_pairings(r1, r2)
    assert result.observation is not None
    assert result.observation.label == "84+_1_1_7_2"
