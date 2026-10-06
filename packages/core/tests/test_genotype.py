"""Genotype caller: model maths, then calls on simulated samples."""

import json
import math
from functools import cache
from pathlib import Path

import numpy as np
import pytest
from scalehd.calibration import HTT_MISEQ, StutterCurve, logit
from scalehd.cli import main
from scalehd.counts import SampleCounts, count_reads
from scalehd.genotype import CallerSettings, Candidate, Flag, NoMoleculesError, call_genotype
from scalehd.genotype.model.call import _same_configuration
from scalehd.genotype.model.kernel import _log_kernel, _log_survival
from scalehd.genotype.model.likelihood import (
    _LONG_ALLELE_SPAN,
    _Model,
    _prior,
    _row_loglik,
    _total_loglik,
    _unpack,
)
from scalehd.genotype.model.table import _Table
from scalehd.simulate import SimAllele, SimulationSpec, simulate
from scalehd.structure import AlleleStructure
from scipy.special import logsumexp

# A sample's stutter moved this far, on each of the six scales of a StutterCurve, from
# the priors' curve.
DRIFT_36 = (0.0336, 0.3961, -1.4568, 2.0608, -1.3735, 0.2539)
DRIFT_131 = (0.029, -0.4377, 0.3564, 0.4743, -0.4999, 1.3539)


def drifted(offsets: tuple[float, ...]) -> StutterCurve:
    """HTT_MISEQ moved by ``offsets``. (N-1)/N and (N+1)/N stay at most 0.95, so N is still
    its allele's tallest peak, and the step and tail ratios stay within 0.001 to 0.97."""
    columns = (
        HTT_MISEQ.log_contraction,
        HTT_MISEQ.logit_contraction_step,
        HTT_MISEQ.logit_contraction_tail,
        HTT_MISEQ.log_expansion,
        HTT_MISEQ.logit_expansion_step,
        HTT_MISEQ.logit_expansion_tail,
    )
    moved = []
    for k, (column, offset) in enumerate(zip(columns, offsets, strict=True)):
        low, high = (-math.inf, math.log(0.95)) if k in (0, 3) else (logit(1e-3), logit(0.97))
        moved.append(tuple(float(min(max(v + offset, low), high)) for v in column))
    return StutterCurve(HTT_MISEQ.cag_lengths, *moved, spread=HTT_MISEQ.spread)


@cache
def sample(
    *labels: str,
    pairs: int = 3000,
    seed: int = 5,
    abundance: float = 1.0,
    drift: tuple[float, ...] | None = None,
) -> SampleCounts:
    """Counts for simulated alleles. The last allele's abundance can be lowered, and the
    stutter expected by the model (for CAG length X) can be modified with ``drift``."""
    alleles = [SimAllele(AlleleStructure.from_label(label)) for label in labels]
    if abundance != 1.0:
        alleles[-1] = SimAllele(alleles[-1].structure, abundance)
    stutter = HTT_MISEQ if drift is None else drifted(drift)
    simulated = simulate(SimulationSpec(tuple(alleles), pairs=pairs, stutter=stutter, seed=seed))
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


def _candidate(label: str) -> Candidate:
    return Candidate(AlleleStructure.from_label(label))


def test_lengths_settled_by_their_own_peak_may_differ_others_must_match() -> None:
    called = _Model.of(_candidate("44_1_1_7_2"), _candidate("46_1_1_7_2"))
    other = _Model.of(_candidate("43_1_1_7_2"), _candidate("45_1_1_7_2"))
    assert _same_configuration(other, called, (2, 2))
    assert not _same_configuration(other, called, (2, 0))
    assert not _same_configuration(other, called, (0, 0))


def test_close_alleles_are_only_as_sure_as_their_own_lengths() -> None:
    """Two expanded alleles two CAG apart, PCR stutter not as expected in model.
    assert appropriate caution in genotyping and what the calls are """
    call = call_genotype(sample("43_1_1_7_2", "45_1_1_7_2", pairs=1519, seed=36, drift=DRIFT_36))
    assert call.posterior < 0.99
    genotypes = [call.label] + [label for label, _ in call.alternatives]
    assert "43_1_1_7_2/45_1_1_7_2" in genotypes


def test_exact_length_is_not_pulled_by_the_edge_of_its_peak() -> None:
    """Stutter different to model again. For long alleles test that effect on genotype.
    Only the reads near the peak decide N, and leaving the rest out shouldn't favour e.g. N+1"""
    call = call_genotype(sample("18_1_1_7_2", "57_1_1_7_2", pairs=2545, seed=131, drift=DRIFT_131))
    assert call.label == "18_1_1_7_2/57_1_1_7_2"


@pytest.mark.parametrize(
    ("labels", "close"),
    [
        # Long alleles stutter more, slippage could be a fair % of another peak candidate
        (("44_1_1_7_2", "46_1_1_7_2"), True),
        # Short alleles stutter little, so two apart their peaks barely overlap.
        (("17_1_1_7_2", "19_1_1_7_2"), False),
        (("19_1_1_7_2", "44_1_1_7_2"), False),
        # Different CCG keeps each allele's molecules apart.
        (("17_1_1_7_2", "19_1_1_10_2"), False),
    ],
)
def test_close_alleles_are_flagged_when_their_stutter_overlaps(
    labels: tuple[str, ...], close: bool
) -> None:
    flags = call_genotype(sample(*labels)).flags
    assert (Flag.CLOSE_ALLELES in flags) is close
    assert Flag.NEIGHBOURING not in flags


def test_neighbouring_alleles_are_not_also_close() -> None:
    flags = call_genotype(sample("17_1_1_7_2", "18_1_1_7_2")).flags
    assert Flag.NEIGHBOURING in flags
    assert Flag.CLOSE_ALLELES not in flags


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
