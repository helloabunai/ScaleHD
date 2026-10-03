"""The worker workload i.e. FASTQ (or a simulation) in, counts and call out."""

import json
from pathlib import Path

import pytest
from scalehd.simulate import SimAllele, SimulationSpec, simulate
from scalehd.structure import AlleleStructure
from scalehd_server.schemas import GenotypeMethod, JobSettings
from scalehd_server.worker import SampleTask, run_sample, run_task

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
    assert call.label == "17_1_1_7_2/43_1_1_7_2"
    assert json.loads((out / "call.json").read_text()) == call.to_dict()
    assert json.loads((out / "counts.json").read_text())["schema"] == "scalehd.counts/2"


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


def _task(tmp_path: Path, truth: str) -> SampleTask:
    return SampleTask(
        sample_id=7,
        name="expanded",
        folder=tmp_path / "expanded",
        settings=JobSettings(method=MODEL),
        simulation={"alleles": ["17_1_1_7_2", "43_1_1_7_2"], "pairs": 2000, "seed": 1},
        truth=truth,
    )


def test_run_task_simulates_calls_and_checks_the_truth(tmp_path: Path) -> None:
    result = run_task(_task(tmp_path, "17_1_1_7_2/43_1_1_7_2"))
    assert result.sample_id == 7
    assert result.genotype == "17_1_1_7_2/43_1_1_7_2"
    assert result.matches_truth is True
    assert result.confidence is not None
    assert result.confidence > 0
    assert all(isinstance(flag, str) for flag in result.flags)
    folder = tmp_path / "expanded"
    assert (folder / "input" / "expanded_R1.fastq.gz").exists()
    assert (folder / "counts.json").exists()
    assert json.loads((folder / "call.json").read_text()) == result.call


def test_run_task_reports_a_call_that_misses_the_truth(tmp_path: Path) -> None:
    assert run_task(_task(tmp_path, "17_1_1_7_2/44_1_1_7_2")).matches_truth is False


def test_run_task_needs_inputs_or_a_recipe(tmp_path: Path) -> None:
    task = SampleTask(
        sample_id=1, name="empty", folder=tmp_path, settings=JobSettings(method=MODEL)
    )
    with pytest.raises(ValueError, match="has no input files"):
        run_task(task)
