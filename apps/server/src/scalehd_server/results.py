"""Individual sample results!! built from its saved counts and genotype call"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any, Literal

from scalehd.counts import SampleCounts
from scalehd.structure import FieldStatus

from .models import Job, Sample
from .schemas import CagBar, CagChart, CcgBar, Cell, JobTag, Reads, SampleDetail, SampleOut
from .workspace import sample_folders

FileKind = Literal["call", "counts", "r1", "r2"]
FILE_KINDS: tuple[FileKind, ...] = ("call", "counts", "r1", "r2")

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


def _cag_chart(counts: SampleCounts, rest: _Rest, alleles: list[dict[str, Any]]) -> CagChart:
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
        files=files,
        previous_id=ids[at - 1] if at > 0 else None,
        next_id=ids[at + 1] if at + 1 < len(ids) else None,
    )
