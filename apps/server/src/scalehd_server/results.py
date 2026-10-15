"""Individual sample results!! built from its saved counts and genotype call"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any, Literal

from scalehd.counts import SampleCounts
from scalehd.structure import FieldStatus

from . import legacy
from .models import Job, Sample
from .schemas import (
    CagBar,
    CagChart,
    CcgBar,
    Cell,
    GenotypeMethod,
    JobSettings,
    JobTag,
    LegacyCagChart,
    LegacyCcgChart,
    LegacySample,
    Reads,
    SampleDetail,
    SampleOut,
)
from .workspace import sample_folders

# LEGACY: 1.x's result.json, its PDF of graphs, and its log.
FileKind = Literal["call", "counts", "r1", "r2", "legacy-result", "legacy-graphs", "legacy-log"]
FILE_KINDS: tuple[FileKind, ...] = (
    "call",
    "counts",
    "r1",
    "r2",
    "legacy-result",
    "legacy-graphs",
    "legacy-log",
)
# 1.x failed in these after genotyping: its genotype stands.
_AFTER_GENOTYPING = ("SNPCalling", "Report/Graph")

# Everything about a structure but its CAG count: CAACAG, CCGCCA, CCG, CCT.
_Rest = tuple[int, int, int, int]


def cag_charts(counts: SampleCounts, call: dict[str, Any] | None) -> list[CagChart]:
    """A CAG distribution for each structure of the called alleles.

    Two alleles that differ only in CAG (i.e. same intervening/CCG) share one chart.
    """
    if call is None:
        return []
    by_rest: dict[_Rest, list[dict[str, Any]]] = {}
    for allele in call["alleles"]:
        rest = (allele["caacag"], allele["ccgcca"], allele["ccg"], allele["cct"])
        by_rest.setdefault(rest, []).append(allele)
    return [_cag_chart(counts, rest, alleles) for rest, alleles in by_rest.items()]


def cag_tallies(counts: SampleCounts, rest: _Rest) -> tuple[Counter[int], Counter[int]]:
    exact: Counter[int] = Counter()
    for structure, n in counts.complete.items():
        if structure.counts[1:] == rest:
            exact[structure.cag] += n
    lower_bounds: Counter[int] = Counter()
    for observation, n in counts.partial.items():
        cag_status, *rest_status = observation.status
        if not (
            all(status is FieldStatus.EXACT for status in rest_status)
            and tuple(observation.counts[1:]) == rest
        ):
            continue
        if cag_status is FieldStatus.LOWER_BOUND:
            lower_bounds[observation.counts[0]] += n
        elif cag_status is FieldStatus.UNCONFIRMED:
            # Read to the tract's end, just not confirmed: almost always that length.
            exact[observation.counts[0]] += n
    return exact, lower_bounds


def _cag_chart(counts: SampleCounts, rest: _Rest, alleles: list[dict[str, Any]]) -> CagChart:
    exact, lower_bounds = cag_tallies(counts, rest)
    called = sorted({allele["cag"] for allele in alleles})
    cags = set(exact) | set(lower_bounds) | set(called)
    return CagChart(
        caacag=rest[0],
        ccgcca=rest[1],
        ccg=rest[2],
        cct=rest[3],
        alleles=sorted({allele["structure"] for allele in alleles}),
        called=called,
        bars=[
            CagBar(cag=cag, molecules=exact[cag], lower_bound=lower_bounds[cag])
            for cag in range(min(cags), max(cags) + 1)
        ],
    )


def ccg_distribution(counts: SampleCounts) -> list[CcgBar]:
    """Complete molecules at each CCG length, across every structure."""
    by_ccg: Counter[int] = Counter()
    for structure, n in counts.complete.items():
        by_ccg[structure.ccg] += n
    if not by_ccg:
        return []
    return [CcgBar(ccg=ccg, molecules=by_ccg[ccg]) for ccg in range(min(by_ccg), max(by_ccg) + 1)]


def cag_ccg_cells(counts: SampleCounts) -> list[Cell]:
    """Complete molecules by CAG and CCG length, non-empty cells only."""
    cells: Counter[tuple[int, int]] = Counter()
    for structure, n in counts.complete.items():
        cells[structure.cag, structure.ccg] += n
    return [Cell(cag=cag, ccg=ccg, molecules=n) for (cag, ccg), n in sorted(cells.items())]


def reads(counts: SampleCounts) -> Reads:
    return Reads(
        molecules=counts.molecules,
        complete=sum(counts.complete.values()),
        partial=sum(counts.partial.values()),
        dropped=counts.dropped,
        unusable=counts.unusable,
        read_outcomes=dict(counts.read_outcomes),
        discordant=dict(counts.discordant),
    )


def sample_file(job: Job, sample: Sample, kind: FileKind) -> Path | None:
    """locate sample file for future download feature"""
    folder = sample_folders(job)[sample.id]
    if kind == "call":
        return folder / "call.json"
    if kind == "counts":
        return folder / "counts.json"
    if kind == "legacy-result":
        return folder / legacy.FOLDER / legacy.RESULT
    if kind == "legacy-log":
        return folder / legacy.LOG
    if kind == "legacy-graphs":
        result = legacy.read_result(folder)
        graphs = legacy.graphs(result) if result else None
        # Only ever a file in the sample's own folder.
        return graphs if graphs and graphs.resolve().is_relative_to(folder.resolve()) else None
    read = 1 if kind == "r1" else 2
    if sample.simulation is not None:
        return folder / "input" / f"{folder.name}_R{read}.fastq.gz"
    path = sample.r1 if read == 1 else sample.r2
    return Path(path) if path else None


def sample_detail(job: Job, sample: Sample) -> SampleDetail:
    """sample result page data"""
    folder = sample_folders(job)[sample.id] if job.output_dir else None
    counts_file = folder / "counts.json" if folder else None
    counts = SampleCounts.read_json(counts_file) if counts_file and counts_file.is_file() else None
    files = (
        [kind for kind in FILE_KINDS if (path := sample_file(job, sample, kind)) and path.is_file()]
        if folder
        else []
    )
    ids = [other.id for other in job.samples]
    at = ids.index(sample.id)
    is_legacy = JobSettings.model_validate(job.settings).method is GenotypeMethod.LEGACY
    return SampleDetail(
        sample=SampleOut.model_validate(sample),
        job_id=job.id,
        job_name=job.name,
        demo=job.demo,
        tags=[JobTag.model_validate(tag) for tag in job.tags],
        folder=str(folder) if folder else None,
        call=sample.call,
        cag_charts=cag_charts(counts, sample.call) if counts else [],
        ccg=ccg_distribution(counts) if counts else [],
        cells=cag_ccg_cells(counts) if counts else [],
        reads=reads(counts) if counts else None,
        legacy=legacy_sample(folder) if folder and is_legacy else None,
        files=files,
        previous_id=ids[at - 1] if at > 0 else None,
        next_id=ids[at + 1] if at + 1 < len(ids) else None,
    )


def legacy_sample(folder: Path) -> LegacySample | None:
    result = legacy.read_result(folder)
    if result is None:
        return None
    genotyped = "genotyping" in result.get("stages", []) and (
        result["status"] == "complete" or result.get("failed_stage") in _AFTER_GENOTYPING
    )
    alleles = [
        (result.get("sample") or {}).get(key) for key in ("primaryallele", "secondaryallele")
    ]
    called = genotyped and all(allele and allele.get("reflabel") for allele in alleles)
    return LegacySample(
        status=result["status"],
        failed_stage=result.get("failed_stage"),
        error=result.get("error"),
        stages=result.get("stages", []),
        seconds=result.get("seconds"),
        threads=result.get("threads"),
        versions=result.get("versions"),
        call=legacy.call(result) if called else None,
        cag_charts=_legacy_cag_charts(alleles) if called else [],
        ccg_charts=_legacy_ccg_charts(alleles) if called else [],
    )


def _legacy_cag_charts(alleles: list[Any]) -> list[LegacyCagChart]:
    charts: list[LegacyCagChart] = []
    for allele in alleles:
        structure = allele["reflabel"]
        cag = int(structure.split("_")[0])
        reads = legacy.cag_reads(allele)
        same = next((c for c in charts if c.ccg == allele["ccg"] and c.reads == reads), None)
        if same is None:
            charts.append(
                LegacyCagChart(alleles=[structure], called=[cag], ccg=allele["ccg"], reads=reads)
            )
        else:
            same.alleles = sorted({*same.alleles, structure})
            same.called = sorted({*same.called, cag})
    return charts


def _legacy_ccg_charts(alleles: list[Any]) -> list[LegacyCcgChart]:
    charts: list[LegacyCcgChart] = []
    for allele in alleles:
        structure, reads = allele["reflabel"], allele["rvarray"]
        same = next((c for c in charts if c.reads == reads), None)
        if same is None:
            charts.append(LegacyCcgChart(alleles=[structure], reads=reads))
        elif structure not in same.alleles:
            same.alleles.append(structure)
    return charts
