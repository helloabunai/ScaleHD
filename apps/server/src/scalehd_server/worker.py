"""What a worker process does with one sample.

A task calls this w/ a plain description of the sample, and a plain
result goes back to the server process. db writing happens there.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from scalehd.counts import count_fastq
from scalehd.genotype import GenotypeCall, call_genotype
from scalehd.simulate import SimAllele, SimulationSpec, call_matches, simulate
from scalehd.structure import AlleleStructure

from . import legacy
from .config import ServerSettings
from .progress import ProgressReporter, SampleStage
from .schemas import GenotypeMethod, JobSettings, LegacySettings


@dataclass(frozen=True)
class SampleTask:
    sample_id: int
    name: str
    folder: Path
    settings: JobSettings
    r1: Path | None = None
    r2: Path | None = None
    # {"alleles": [labels], "pairs": n, "seed": s} for a simulated sample.
    simulation: dict[str, Any] | None = None
    # The true genotype label, for a simulated sample.
    truth: str | None = None
    threads: int = 1
    # Runs ScaleHD 1.x (ServerSettings.legacy_runner)
    legacy_runner: Path | None = None


@dataclass(frozen=True)
class SampleResult:
    sample_id: int
    call: dict[str, Any] | None
    genotype: str | None
    # The call's Phred-scaled confidence (GenotypeCall.quality).
    confidence: float | None
    flags: list[str]
    matches_truth: bool | None


def run_sample(
    r1: Path,
    r2: Path | None,
    out_dir: Path,
    settings: JobSettings,
    progress: ProgressReporter | None = None,
) -> GenotypeCall | None:
    """Count one sample's molecules and call its genotype.

    Writes ``counts.json`` and, when calling, ``call.json`` to out_dir. Returns the
    call, or None for a count-only job.
    """
    if settings.method is not GenotypeMethod.MODEL:
        raise ValueError(f"{settings.method} genotyping doesn't count reads")
    out_dir.mkdir(parents=True, exist_ok=True)
    if progress is not None:
        progress.stage(SampleStage.COUNTING, 0.0)
    counted = progress.counted if progress is not None else None
    counts = count_fastq(r1, r2, policy=settings.discordant, progress=counted)
    counts.write_json(out_dir / "counts.json")
    if not settings.call:
        return None
    if progress is not None:
        progress.stage(SampleStage.CALLING)
    call = call_genotype(counts, settings.caller_settings())
    (out_dir / "call.json").write_text(json.dumps(call.to_dict(), indent=2) + "\n")
    return call


def run_task(task: SampleTask) -> SampleResult:
    """Simulate the sample if it is simulated, then count and call it."""
    progress = ProgressReporter(task.folder)
    r1, r2 = task.r1, task.r2
    if task.simulation is not None:
        progress.stage(SampleStage.SIMULATING)
        r1, r2 = _simulate(task.folder, task.simulation)
    if r1 is None:
        raise ValueError(f"sample {task.name} has no input files")
    if task.settings.method is GenotypeMethod.LEGACY:
        return _run_legacy(task, r1, r2, progress)
    call = run_sample(r1, r2, task.folder, task.settings, progress)
    if call is None:
        return SampleResult(task.sample_id, None, None, None, [], None)
    matches = None
    if task.truth is not None:
        truth = [AlleleStructure.from_label(label) for label in task.truth.split("/")]
        matches = call_matches([allele.allele for allele in call.alleles], truth)
    return SampleResult(
        sample_id=task.sample_id,
        call=call.to_dict(),
        genotype=call.label,
        confidence=call.quality,
        flags=[str(flag) for flag in call.flags],
        matches_truth=matches,
    )


def available_methods(config: ServerSettings) -> list[GenotypeMethod]:
    """The genotyping methods this server can run. Legacy needs ScaleHD 1.x, which only
    the Docker image has."""
    return [
        method
        for method in GenotypeMethod
        if method is not GenotypeMethod.LEGACY or config.legacy_available
    ]


def _run_legacy(
    task: SampleTask, r1: Path, r2: Path | None, progress: ProgressReporter
) -> SampleResult:
    if task.legacy_runner is None:
        raise legacy.LegacyError("this server can't run legacy ScaleHD")
    result = legacy.run_sample(
        task.name,
        task.folder,
        r1,
        r2,
        task.settings.legacy or LegacySettings(),
        task.threads,
        task.legacy_runner,
        progress,
    )
    call = legacy.call(result)
    return SampleResult(
        sample_id=task.sample_id,
        call=call,
        genotype=call["genotype"],
        confidence=legacy.confidence(call),
        flags=call["flags"],
        matches_truth=legacy.matches(call, task.truth) if task.truth is not None else None,
    )


def _simulate(folder: Path, recipe: dict[str, Any]) -> tuple[Path, Path]:
    alleles = tuple(SimAllele(AlleleStructure.from_label(label)) for label in recipe["alleles"])
    spec = SimulationSpec(alleles, pairs=recipe["pairs"], seed=recipe["seed"])
    r1, r2, _ = simulate(spec).write(folder / "input", folder.name)
    return r1, r2
