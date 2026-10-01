from scalehd.pairs import DiscordancePolicy, join_mates
from scalehd.structure import AlleleStructure, FieldStatus, Observation

E, L, U = FieldStatus.EXACT, FieldStatus.LOWER_BOUND, FieldStatus.UNOBSERVED


def exact(label: str) -> Observation:
    return Observation.exact(AlleleStructure.from_label(label))


def test_agreeing_mates() -> None:
    result = join_mates(exact("43_1_1_7_2"), exact("43_1_1_7_2"))
    assert result.observation == exact("43_1_1_7_2")
    assert not result.is_discordant


def test_disagreeing_mates_dropped_by_default() -> None:
    result = join_mates(exact("43_1_1_7_2"), exact("41_2_1_7_2"))
    assert result.observation is None
    assert result.discordant == (True, True, False, False, False)


def test_prefer_policy_takes_preferred_mate_per_field() -> None:
    result = join_mates(exact("43_1_1_7_2"), exact("43_1_1_8_2"), policy=DiscordancePolicy.PREFER)
    assert result.observation is not None
    assert result.observation.label == "43_1_1_8_2"  # CCG comes from R2


def test_long_allele_completed_from_both_mates() -> None:
    r1 = Observation((75, 1, 1, 2, 0), (E, E, E, L, U))
    r2 = Observation((60, 1, 1, 7, 2), (L, E, E, E, E))
    result = join_mates(r1, r2)
    assert result.observation is not None
    assert result.observation.label == "75_1_1_7_2"
    assert not result.is_discordant


def test_censored_when_neither_mate_spans_cag() -> None:
    r1 = Observation((84, 0, 0, 0, 0), (L, U, U, U, U))
    r2 = Observation((70, 1, 1, 7, 2), (L, E, E, E, E))
    result = join_mates(r1, r2)
    assert result.observation is not None
    assert result.observation.label == "84+_1_1_7_2"


def test_lower_bound_above_exact_is_discordant() -> None:
    r1 = Observation((40, 1, 1, 7, 2), (E, E, E, E, E))
    r2 = Observation((45, 1, 1, 7, 2), (L, E, E, E, E))
    assert join_mates(r1, r2).observation is None


def test_single_mate() -> None:
    assert join_mates(exact("20_1_1_7_2"), None).observation == exact("20_1_1_7_2")
    assert join_mates(None, None).observation is None
