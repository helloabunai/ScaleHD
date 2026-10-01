"""End to end: simulated FASTQ through parsing and read pair joining."""

import json
from pathlib import Path

import pytest
from scalehd.counts import SampleCounts, count_fastq, count_reads
from scalehd.pairs import join_mates
from scalehd.parse import RepeatParser
from scalehd.seqio import read_pairs, reverse_complement
from scalehd.simulate import SimAllele, SimulationSpec, StutterModel, simulate
from scalehd.structure import AlleleStructure, FieldStatus


def spec(*labels: str, pairs: int = 3000, seed: int = 7, **kwargs: object) -> SimulationSpec:
    alleles = tuple(SimAllele(AlleleStructure.from_label(label)) for label in labels)
    return SimulationSpec(alleles, pairs=pairs, seed=seed, **kwargs)  # type: ignore[arg-type]


def test_simulation_is_deterministic() -> None:
    a = simulate(spec("17_1_1_7_2", "43_1_1_7_2", pairs=200))
    b = simulate(spec("17_1_1_7_2", "43_1_1_7_2", pairs=200))
    assert a.r1 == b.r1
    assert a.r2 == b.r2
    assert a.molecules == b.molecules


def test_reads_have_requested_shape() -> None:
    sample = simulate(spec("20_1_1_7_2", pairs=50))
    assert len(sample.r1) == len(sample.r2) == 50
    assert all(len(r.sequence) == 300 == len(r.quality) for r in sample.r1 + sample.r2)
    assert sum(sample.molecules.values()) == 50


def test_stutter_probabilities_scale_with_length() -> None:
    model = StutterModel()
    short, long = model.cag_shift_probabilities(20), model.cag_shift_probabilities(60)
    assert short.sum() == pytest.approx(1)
    assert long[1] > short[1]


def test_every_kept_molecule_matches_truth() -> None:
    """With pair agreement, sequencing errors should not produce wrong molecules."""
    sample = simulate(spec("17_1_1_7_2", "43_1_1_7_2", "42_0_1_7_2"))
    parser = RepeatParser()
    kept = wrong = 0
    for a, b in zip(sample.r1, sample.r2, strict=True):
        x = parser.parse(a.sequence)
        y = parser.parse(reverse_complement(b.sequence))
        joined = join_mates(
            x.observation if x.usable else None, y.observation if y.usable else None
        )
        if joined.observation is None:
            continue
        kept += 1
        wrong += joined.observation.label != a.name.rsplit(":", 1)[-1]
    assert kept > 0.9 * len(sample.r1)
    assert wrong <= 0.001 * kept


@pytest.mark.parametrize(
    "labels",
    [
        ("17_1_1_7_2", "43_1_1_7_2"),
        ("21_1_1_7_2",),
        ("42_0_1_7_2", "19_2_1_10_2"),
        ("20_1_1_7_2", "75_1_1_7_2"),
    ],
)
def test_true_alleles_are_the_top_structures(labels: tuple[str, ...]) -> None:
    sample = simulate(spec(*labels))
    counts = count_reads(
        (a.sequence, b.sequence) for a, b in zip(sample.r1, sample.r2, strict=True)
    )
    top = {s.label for s, _ in counts.top(len(labels))}
    assert top == set(labels)


def test_very_long_allele_is_censored_not_miscalled() -> None:
    sample = simulate(spec("20_1_1_7_2", "95_1_1_7_2"))
    counts = count_reads(
        (a.sequence, b.sequence) for a, b in zip(sample.r1, sample.r2, strict=True)
    )
    assert counts.top(1)[0][0].label == "20_1_1_7_2"
    assert not any(s.cag > 70 for s in counts.complete)
    observation, _ = counts.partial.most_common(1)[0]
    assert observation.status[0] is FieldStatus.LOWER_BOUND
    assert observation.label.endswith("+_1_1_7_2")


def test_off_target_pairs_are_unusable() -> None:
    sample = simulate(spec("20_1_1_7_2", pairs=500, off_target=0.2))
    counts = count_reads(
        (a.sequence, b.sequence) for a, b in zip(sample.r1, sample.r2, strict=True)
    )
    assert counts.unusable == sample.off_target_pairs > 0


def test_files_round_trip(tmp_path: Path) -> None:
    sample = simulate(spec("17_1_1_7_2", "43_1_1_7_2", pairs=300))
    r1, r2, truth = sample.write(tmp_path, "s1")
    assert [x.sequence for x, _ in read_pairs(r1, r2)] == [x.sequence for x in sample.r1]
    assert json.loads(truth.read_text())["alleles"][1]["structure"] == "43_1_1_7_2"

    counts = count_fastq(r1, r2)
    out = tmp_path / "s1.counts.json"
    counts.write_json(out)
    assert SampleCounts.read_json(out) == counts
    matrix = counts.cag_ccg_matrix()
    assert matrix.shape == (20, 200)
    assert matrix.sum() == sum(counts.complete.values())
    assert matrix[6, 16] == counts.complete[AlleleStructure(17)]
