"""Genotype caller: model maths, then calls on simulated samples."""

import json
from collections import Counter
from functools import cache
from pathlib import Path

import numpy as np
import pytest
from scalehd.calibration import HTT_MISEQ
from scalehd.cli import main
from scalehd.counts import SampleCounts, count_reads
from scalehd.genotype import (
    Flag,
    NoMoleculesError,
    _log_kernel,
    _log_survival,
    call_genotype,
    coarsen,
)
from scalehd.simulate import SimAllele, SimulationSpec, simulate
from scalehd.structure import AlleleStructure, FieldStatus, Observation


@cache
def sample(*labels: str, pairs: int = 3000, seed: int = 5, abundance: float = 1.0) -> SampleCounts:
    """Counts for simulated alleles; the last allele's abundance can be lowered."""
    alleles = [SimAllele(AlleleStructure.from_label(label)) for label in labels]
    if abundance != 1.0:
        alleles[-1] = SimAllele(alleles[-1].structure, abundance)
    simulated = simulate(SimulationSpec(tuple(alleles), pairs=pairs, seed=seed))
    return count_reads(
        (a.sequence, b.sequence) for a, b in zip(simulated.r1, simulated.r2, strict=True)
    )


@pytest.mark.parametrize("n", [1, 2, 3, 17, 44, 80])
def test_kernel_normalises_and_survival_matches_it(n: int) -> None:
    stutter = HTT_MISEQ.at(n)
    values = np.arange(1, n + 400)
    probabilities = np.exp(_log_kernel(values, n, stutter))
    assert probabilities.sum() == pytest.approx(1.0)
    bounds = np.array([1, 2, max(n - 3, 1), n, n + 1, n + 5])
    expected = [probabilities[values >= b].sum() for b in bounds]
    assert np.exp(_log_survival(bounds, n, stutter)) == pytest.approx(expected, abs=1e-9)


def test_coarsen_keeps_every_molecule() -> None:
    counts = SampleCounts(molecules=6)
    counts.complete[AlleleStructure(20)] = 3
    counts.complete[AlleleStructure(80)] = 2
    bounded = (FieldStatus.LOWER_BOUND,) + (FieldStatus.EXACT,) * 4
    counts.partial[Observation((83, 1, 1, 7, 2), bounded)] = 1
    out = coarsen(counts, 77)
    assert out.complete == Counter({AlleleStructure(20): 3})
    assert out.partial == Counter({Observation((78, 1, 1, 7, 2), bounded): 3})


@pytest.mark.parametrize(
    ("labels", "expected"),
    [
        (("17_1_1_7_2", "43_1_1_7_2"), "17_1_1_7_2/43_1_1_7_2"),
        (("21_1_1_7_2",), "21_1_1_7_2/21_1_1_7_2"),
        (("17_1_1_7_2", "18_1_1_7_2"), "17_1_1_7_2/18_1_1_7_2"),
        (("16_1_1_7_2", "17_1_1_7_2"), "16_1_1_7_2/17_1_1_7_2"),
        (("17_1_1_7_2", "17_1_1_10_2"), "17_1_1_7_2/17_1_1_10_2"),
        (("42_0_1_7_2", "19_2_1_10_2"), "19_2_1_10_2/42_0_1_7_2"),
        (("19_1_1_7_2", "44_1_1_7_2"), "19_1_1_7_2/44_1_1_7_2"),
    ],
)
def test_calls_simulated_genotypes(labels: tuple[str, ...], expected: str) -> None:
    call = call_genotype(sample(*labels))
    assert call.label == expected
    assert call.posterior > 0.9


def test_flags_describe_the_genotype() -> None:
    assert Flag.HOMOZYGOUS in call_genotype(sample("21_1_1_7_2")).flags
    assert Flag.NEIGHBOURING in call_genotype(sample("17_1_1_7_2", "18_1_1_7_2")).flags
    atypical = call_genotype(sample("42_0_1_7_2", "19_2_1_10_2"))
    assert Flag.ATYPICAL in atypical.flags
    assert Flag.HOMOZYGOUS not in atypical.flags


def test_allele_beyond_read_length_is_a_lower_bound() -> None:
    call = call_genotype(sample("20_1_1_7_2", "95_1_1_7_2"))
    short, long = call.alleles
    assert short.label == "20_1_1_7_2"
    assert long.allele.beyond_read_length
    assert long.label.endswith("+_1_1_7_2")
    assert long.allele.structure.cag <= 95
    assert Flag.BEYOND_READ_LENGTH in call.flags
    assert long.cag_estimate is not None
    _, low, high = long.cag_estimate
    assert low <= 95 <= high


def test_low_depth_is_less_certain() -> None:
    deep = call_genotype(sample("17_1_1_7_2", "18_1_1_7_2"))
    shallow = call_genotype(sample("17_1_1_7_2", "18_1_1_7_2", pairs=150))
    assert shallow.quality < deep.quality
    assert Flag.LOW_DEPTH in shallow.flags


def test_third_allele_is_reported_as_unexplained() -> None:
    call = call_genotype(sample("17_1_1_7_2", "43_1_1_7_2", "25_1_1_7_2", abundance=0.15))
    assert Flag.UNEXPLAINED_PEAK in call.flags
    assert any(s.startswith("25_") for s, _ in call.unexplained)


def test_call_serialises_to_json() -> None:
    call = call_genotype(sample("17_1_1_7_2", "43_1_1_7_2"))
    data = json.loads(json.dumps(call.to_dict()))
    assert data["genotype"] == "17_1_1_7_2/43_1_1_7_2"
    assert [a["cag"] for a in data["alleles"]] == [17, 43]
    assert data["alleles"][1]["polyglutamine_length"] == 45


def test_empty_sample_is_an_error() -> None:
    with pytest.raises(NoMoleculesError):
        call_genotype(SampleCounts())


def test_cli_call(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    counts_path, call_path = tmp_path / "s.counts.json", tmp_path / "s.call.json"
    sample("17_1_1_7_2", "43_1_1_7_2").write_json(counts_path)
    assert main(["call", str(counts_path), "-o", str(call_path)]) == 0
    assert "17_1_1_7_2/43_1_1_7_2" in capsys.readouterr().out
    assert json.loads(call_path.read_text())["schema"] == "scalehd.call/1"
