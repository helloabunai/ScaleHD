"""Genotype caller: model maths, then calls on simulated samples."""

import json
import math
from functools import cache
from pathlib import Path

import numpy as np
import pytest
from scalehd.calibration import HTT_MISEQ
from scalehd.cli import main
from scalehd.counts import SampleCounts, count_reads
from scalehd.genotype import (
    _LONG_ALLELE_SPAN,
    CallerSettings,
    Candidate,
    Flag,
    NoMoleculesError,
    _log_kernel,
    _log_survival,
    _Model,
    _prior,
    _row_loglik,
    _Table,
    _total_loglik,
    _unpack,
    call_genotype,
)
from scalehd.simulate import SimAllele, SimulationSpec, simulate
from scalehd.structure import AlleleStructure
from scipy.special import logsumexp


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
        (("19_1_0_7_2", "40_1_1_7_2"), "19_1_0_7_2/40_1_1_7_2"),
        (("17_1_1_7_2", "42_1_2_7_2"), "17_1_1_7_2/42_1_2_7_2"),
        # Long alleles, whose contractions reach further below N than short ones'.
        (("20_1_1_7_2", "55_1_1_7_2"), "20_1_1_7_2/55_1_1_7_2"),
        (("20_1_1_7_2", "60_1_1_7_2"), "20_1_1_7_2/60_1_1_7_2"),
        (("20_1_1_7_2", "66_1_1_7_2"), "20_1_1_7_2/66_1_1_7_2"),
        (("20_1_1_7_2", "70_1_1_7_2"), "20_1_1_7_2/70_1_1_7_2"),
        # Near read length (300 bases): many molecules' CAG ends are seen but unconfirmed.
        (("20_1_1_7_2", "76_1_1_7_2"), "20_1_1_7_2/76_1_1_7_2"),
        (("20_1_1_7_2", "80_1_1_7_2"), "20_1_1_7_2/80_1_1_7_2"),
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
    # No read spans the tract, so the reads set no upper limit and there is no estimate.
    assert long.cag_estimate is None


def test_beyond_read_length_allele_has_one_length_for_all_its_molecules() -> None:
    """The average over N of the whole sample's likelihood, not of each molecule's."""
    counts = sample("20_1_1_7_2", "95_1_1_7_2")
    settings = CallerSettings()
    table = _Table(
        counts, settings.effective_molecules, settings.stutter_window, settings.unconfirmed_error
    )
    short = Candidate(AlleleStructure.from_label("20_1_1_7_2"))
    long = Candidate(AlleleStructure.from_label("83_1_1_7_2"), beyond_read_length=True)
    model = _Model.of(short, long)
    params = _unpack(_prior(model, settings)[0], model)
    each_n = [
        _total_loglik(
            table,
            _row_loglik(
                table, _Model.of(short, Candidate(long.structure.with_counts(cag=n))), params
            )[0],
        )
        for n in range(83, 83 + _LONG_ALLELE_SPAN)
    ]
    expected = float(logsumexp(each_n)) - math.log(_LONG_ALLELE_SPAN)
    assert _total_loglik(table, _row_loglik(table, model, params)[0]) == pytest.approx(expected)


@pytest.mark.parametrize("cag", [84, 86])
def test_allele_just_past_read_length_is_a_lower_bound(cag: int) -> None:
    # Reads stop at about 83 CAG, so these can't be told apart from longer ones.
    call = call_genotype(sample("20_1_1_7_2", f"{cag}_1_1_7_2"))
    short, long = call.alleles
    assert short.label == "20_1_1_7_2"
    assert long.allele.beyond_read_length
    assert long.allele.structure.cag <= cag


def test_estimate_for_an_allele_just_past_read_length_covers_it() -> None:
    long = call_genotype(sample("20_1_1_7_2", "90_1_1_7_2")).alleles[1]
    assert long.allele.beyond_read_length
    assert long.cag_estimate is not None
    _, low, high = long.cag_estimate
    assert low <= 90 <= high


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
    assert json.loads(call_path.read_text())["schema"] == "scalehd.call/2"
