"""The job runner's sample step: FASTQ in, counts and call out."""

import json
from pathlib import Path

import pytest
from scalehd.simulate import SimAllele, SimulationSpec, simulate
from scalehd.structure import AlleleStructure
from scalehd_server.runner import run_sample
from scalehd_server.schemas import GenotypeMethod, JobSettings

MODEL = GenotypeMethod.MODEL


def _simulated_pair(directory: Path) -> tuple[Path, Path]:
    alleles = tuple(SimAllele(AlleleStructure.from_label(x)) for x in ("17_1_1_7_2", "43_1_1_7_2"))
    r1, r2, _ = simulate(SimulationSpec(alleles, pairs=2000, seed=1)).write(directory, "s1")
    return r1, r2


def test_run_sample_counts_and_calls(tmp_path: Path) -> None:
    r1, r2 = _simulated_pair(tmp_path / "input")
    out = tmp_path / "out"
    call = run_sample(r1, r2, out, JobSettings(method=MODEL))
    assert call is not None
    assert call["genotype"] == "17_1_1_7_2/43_1_1_7_2"
    assert json.loads((out / "call.json").read_text()) == call
    assert json.loads((out / "counts.json").read_text())["schema"] == "scalehd.counts/1"


def test_count_only_job_skips_the_call(tmp_path: Path) -> None:
    r1, r2 = _simulated_pair(tmp_path / "input")
    out = tmp_path / "out"
    assert run_sample(r1, r2, out, JobSettings(method=MODEL, call=False)) is None
    assert (out / "counts.json").exists()
    assert not (out / "call.json").exists()


def test_job_thresholds_reach_the_caller() -> None:
    caller = JobSettings(min_molecules=50, min_posterior=0.9).caller_settings()
    assert (caller.min_molecules, caller.min_posterior) == (50, 0.9)
    assert caller.max_background == JobSettings().caller_settings().max_background


def test_legacy_method_is_refused_until_extracted(tmp_path: Path) -> None:
    r1 = tmp_path / "s1_R1.fastq"
    with pytest.raises(ValueError, match="legacy genotyping is not available yet"):
        run_sample(r1, None, tmp_path / "out", JobSettings(method=GenotypeMethod.LEGACY))
