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

from .schemas import GenotypeMethod, JobSettings

# Methods jobs can use. SHD 1.x joins once it is extracted from legacy/.
RUNNABLE_METHODS = frozenset({GenotypeMethod.MODEL})


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


@dataclass(frozen=True)
class SampleResult:
    sample_id: int
    call: dict[str, Any] | None
    genotype: str | None
    quality: float | None
    flags: list[str]
    matches_truth: bool | None


def run_sample(
    r1: Path, r2: Path | None, out_dir: Path, settings: JobSettings
) -> GenotypeCall | None:
    """Count one sample's molecules and call its genotype.

    Writes ``counts.json`` and, when calling, ``call.json`` to out_dir. Returns the
    call, or None for a count-only job.
    """
    if settings.method not in RUNNABLE_METHODS:
        raise ValueError(f"{settings.method} genotyping is not available yet")
    out_dir.mkdir(parents=True, exist_ok=True)
    counts = count_fastq(r1, r2, policy=settings.discordant)
    counts.write_json(out_dir / "counts.json")
    if not settings.call:
        return None
    call = call_genotype(counts, settings.caller_settings())
    (out_dir / "call.json").write_text(json.dumps(call.to_dict(), indent=2) + "\n")
    return call


def run_task(task: SampleTask) -> SampleResult:
    """Simulate the sample if it is simulated, then count and call it."""
    r1, r2 = task.r1, task.r2
    if task.simulation is not None:
        r1, r2 = _simulate(task.folder, task.simulation)
    if r1 is None:
        raise ValueError(f"sample {task.name} has no input files")
    call = run_sample(r1, r2, task.folder, task.settings)
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
        quality=call.quality,
        flags=[str(flag) for flag in call.flags],
        matches_truth=matches,
    )


def _simulate(folder: Path, recipe: dict[str, Any]) -> tuple[Path, Path]:
    alleles = tuple(SimAllele(AlleleStructure.from_label(label)) for label in recipe["alleles"])
    spec = SimulationSpec(alleles, pairs=recipe["pairs"], seed=recipe["seed"])
    r1, r2, _ = simulate(spec).write(folder / "input", folder.name)
    return r1, r2
