"""A job's results as CSV tables, a row per sample.

The results export is like ScaleHD 1.x's InstanceReport.csv, for opening in a spreadsheet,
built from the calls stored with each sample.

Validation export is for modifying the statistical model on results processed with real data.
Will be removed after beta or something. Only for development use but you can click it if you
want bruh.
"""

from __future__ import annotations

import csv
import io
import json
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import scalehd
from scalehd.counts import SampleCounts
from scalehd.seqio import read_fastq

from . import __version__
from .models import Job, Sample
from .results import cag_tallies, sample_file
from .schemas import GenotypeMethod, JobSettings
from .workspace import folder_name, sample_folders

# Each allele's columns: the header after "Allele 1 " or "Allele 2 ", and the call's key.
_ALLELE = (
    ("CAG", "cag"),
    ("CAG estimate", "cag_estimate"),
    ("CCG", "ccg"),
    ("share", "fraction"),
    ("molecules", "molecules"),
    ("backward slippage", "backward_slippage"),
    ("somatic mosaicism", "somatic_mosaicism"),
    ("expansion index", "expansion_index"),
    ("contraction index", "contraction_index"),
)
HEADER = [
    "Sample",
    "Validation ID",
    "Status",
    "Genotype",
    "Quality",
    "Posterior",
    "Flags",
    "Molecules",
    *(
        column
        for n in (1, 2)
        for column in (f"Allele {n}", *(f"Allele {n} {name}" for name, _ in _ALLELE))
    ),
    "Truth",
    "Matches truth",
    "Error",
    "R1",
    "R2",
]

# The fitted stutter ratios (scalehd.calibration.Stutter), as headers and keys.
_STUTTER = (
    ("contraction", "contraction"),
    ("contraction step", "contraction_step"),
    ("contraction tail", "contraction_tail"),
    ("expansion", "expansion"),
    ("expansion step", "expansion_step"),
    ("expansion tail", "expansion_tail"),
)
# Molecules around each allele's peak, from N-20 to N+30: the reach of the stutter kernel.
PEAK_WINDOW = (-20, 30)
_READ_OUTCOMES = ("spanning", "truncated", "no_anchor", "nonconforming")
LAB_COLUMNS = ["Lab genotype", "Lab genotype method", "Replicate group", "Sample type", "Lab notes"]
VALIDATION_HEADER = [
    "Validation ID",
    "Batch",
    *LAB_COLUMNS,
    "Status",
    "Genotype",
    "Quality",
    "Posterior",
    "Flags",
    *(column for n in (1, 2, 3) for column in (f"Alternative {n}", f"Alternative {n} posterior")),
    "Background",
    "CCG slippage",
    "Misread",
    "Unexplained peaks",
    *(
        column
        for n in (1, 2)
        for column in (
            f"Allele {n}",
            *(f"Allele {n} {name}" for name, _ in _ALLELE),
            *(f"Allele {n} stutter {name}" for name, _ in _STUTTER),
            f"Allele {n} molecules N{PEAK_WINDOW[0]} to N+{PEAK_WINDOW[1]}",
        )
    ),
    "Molecules",
    "Complete",
    "Partial",
    "Dropped",
    "Unusable",
    *(f"R{read} {outcome.replace('_', ' ')}" for read in (1, 2) for outcome in _READ_OUTCOMES),
    "Discordant",
    "Read length",
    "Truth",
    "Matches truth",
    "Error",
    "ScaleHD version",
    "Server version",
    "Method",
    "Settings",
]


def results_filename(job: Job) -> str:
    label = "demo" if job.demo else folder_name(job.name, "job")
    return f"{job.id}-{label}-results.csv"


def validation_filename(job: Job) -> str:
    # No job name: it might say whose samples these are.
    return f"scalehd-validation-job-{job.id}.csv"


def validation_id(job: Job, sample: Sample) -> str:
    """A sample's ID for sharing: the job's number and the sample's place in it."""
    return f"{job.id}-{job.samples.index(sample) + 1:03d}"


def results_csv(job: Job) -> str:
    return _table(HEADER, (_row(job, sample) for sample in job.samples))


def validation_csv(job: Job) -> str:
    folders = sample_folders(job) if job.output_dir else {}
    versions = _versions(job)
    batches: dict[str, str] = {}
    rows = []
    for sample in job.samples:
        folder = folders.get(sample.id)
        counts_file = folder / "counts.json" if folder else None
        counts = (
            SampleCounts.read_json(counts_file) if counts_file and counts_file.is_file() else None
        )
        rows.append(_validation_row(job, sample, counts, versions, batches))
    return _table(VALIDATION_HEADER, rows)


def _table(header: list[str], rows: Any) -> str:
    out = io.StringIO()
    writer = csv.writer(out, lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return out.getvalue()


def _row(job: Job, sample: Sample) -> list[str]:
    call = sample.call or {}
    row = [
        _text(sample.name),
        validation_id(job, sample),
        str(sample.status),
        *_call_summary(call),
        _number(call.get("molecules"), 0),
    ]
    for allele in _alleles(call):
        row.append(allele.get("structure", ""))
        row.extend(_allele_value(allele, key) for _, key in _ALLELE)
    row += [
        sample.truth or "",
        _yes_no(sample.matches_truth),
        _text(sample.error or ""),
        _text(sample.r1 or ""),
        _text(sample.r2 or ""),
    ]
    return row


def _validation_row(
    job: Job,
    sample: Sample,
    counts: SampleCounts | None,
    versions: dict[str, str],
    batches: dict[str, str],
) -> list[str]:
    call = sample.call or {}
    row = [validation_id(job, sample), _batch(sample, batches), *([""] * len(LAB_COLUMNS))]
    row += [str(sample.status), *_call_summary(call)]
    alternatives = call.get("alternatives", [])[:3]
    for n in range(3):
        alternative = alternatives[n] if n < len(alternatives) else {}
        row += [alternative.get("genotype", ""), _significant(alternative.get("posterior"))]
    row += [
        _significant(call.get("background")),
        _significant(call.get("ccg_slippage")),
        _significant(call.get("misread")),
        "; ".join(f"{p['structure']} ({p['molecules']})" for p in call.get("unexplained", [])),
    ]
    for allele in _alleles(call):
        row.append(allele.get("structure", ""))
        row.extend(_allele_value(allele, key) for _, key in _ALLELE)
        stutter = allele.get("stutter") or {}
        row.extend(_significant(stutter.get(key)) for _, key in _STUTTER)
        row.append(_peak_molecules(counts, allele))
    row += _read_columns(counts)
    row += [
        _read_length(job, sample),
        sample.truth or "",
        _yes_no(sample.matches_truth),
        # Just the kind of error: the message can name a sample or a file.
        (sample.error or "").split(":", 1)[0],
        versions.get("scalehd", ""),
        versions.get("scalehd-server", ""),
    ]
    settings = JobSettings.model_validate(job.settings)
    row += [str(settings.method), _settings(settings)]
    return row


def _call_summary(call: dict[str, Any]) -> list[str]:
    return [
        call.get("genotype", ""),
        _number(call.get("quality"), 1),
        _number(call.get("posterior"), 4),
        "; ".join(call.get("flags", [])),
    ]


def _alleles(call: dict[str, Any]) -> list[dict[str, Any]]:
    """The call's two alleles, or blanks for a sample not called."""
    alleles: list[dict[str, Any]] = call.get("alleles", [])
    return [alleles[n] if n < len(alleles) else {} for n in (0, 1)]


def _allele_value(allele: dict[str, Any], key: str) -> str:
    value = allele.get(key)
    if key == "cag_estimate":
        return f"{value[0]} ({value[1]}-{value[2]})" if value else ""
    if key == "cag" and allele.get("beyond_read_length"):
        # Too long to read: only a lower bound.
        return f"{value}+"
    if key in ("cag", "ccg"):
        return "" if value is None else str(value)
    if key == "molecules":
        return _number(value, 0)
    return _number(value, 3)


def _peak_molecules(counts: SampleCounts | None, allele: dict[str, Any]) -> str:
    """Molecules of the allele's structure at each CAG length around its peak, space
    separated. Blank for an allele beyond read length as its peak was never read."""
    if counts is None or not allele or allele.get("beyond_read_length"):
        return ""
    rest = (allele["caacag"], allele["ccgcca"], allele["ccg"], allele["cct"])
    exact, _ = cag_tallies(counts, rest)
    low, high = PEAK_WINDOW
    return " ".join(str(exact[cag]) for cag in range(allele["cag"] + low, allele["cag"] + high + 1))


def _read_columns(counts: SampleCounts | None) -> list[str]:
    if counts is None:
        return [""] * (5 + 2 * len(_READ_OUTCOMES) + 1)
    outcomes = [
        str(counts.read_outcomes.get(f"r{read}.{outcome}", 0))
        for read in (1, 2)
        for outcome in _READ_OUTCOMES
    ]
    return [
        str(counts.molecules),
        str(sum(counts.complete.values())),
        str(sum(counts.partial.values())),
        str(counts.dropped),
        str(counts.unusable),
        *outcomes,
        "; ".join(f"{name} {n}" for name, n in sorted(counts.discordant.items())),
    ]


def _read_length(job: Job, sample: Sample) -> str:
    if not job.output_dir:
        return ""
    path = sample_file(job, sample, "r1")
    try:
        return str(len(next(read_fastq(path)).sequence)) if path else ""
    except (OSError, ValueError, StopIteration):
        return ""


def _batch(sample: Sample, batches: dict[str, str]) -> str:
    """Origin folder (to be masked)"""
    if not sample.r1:
        return ""
    folder = str(Path(sample.r1).parent)
    return batches.setdefault(folder, f"B{len(batches) + 1}")


def _versions(job: Job) -> dict[str, str]:
    """software version for referencing"""
    running = {"scalehd": scalehd.__version__, "scalehd-server": __version__}
    if not job.output_dir:
        return running
    try:
        recorded = json.loads((Path(job.output_dir) / "job.json").read_text())["versions"]
    except (OSError, ValueError, KeyError, TypeError):
        return running
    return {name: str(recorded.get(name, "")) for name in running}


def _settings(settings: JobSettings) -> str:
    """Used job settings (sharing). A legacy job has ScaleHD 1.x settings. Model settings not
    exposed yet."""
    if settings.method is GenotypeMethod.LEGACY:
        shown = settings.model_dump(mode="json", include={"legacy"})
    else:
        shown = settings.model_dump(mode="json", exclude={"method", "legacy"})
    return "; ".join(
        f"{name}={'default' if value is None else value}" for name, value in _flat(shown)
    )


def _flat(values: dict[str, Any], prefix: str = "") -> Iterator[tuple[str, Any]]:
    """Nested settings as dotted names, e.g. legacy.trim_flags.error_tolerance."""
    for name, value in values.items():
        if isinstance(value, dict):
            yield from _flat(value, f"{prefix}{name}.")
        else:
            yield f"{prefix}{name}", value


def _number(value: float | None, places: int) -> str:
    if value is None:
        return ""
    return f"{value:.{places}f}"


def _significant(value: float | None) -> str:
    """Small fitted values, to 4 significant figures."""
    return "" if value is None else f"{value:.4g}"


def _yes_no(value: bool | None) -> str:
    return {True: "yes", False: "no", None: ""}[value]


def _text(value: str) -> str:
    """Free text, kept from being read as a formula by a spreadsheet."""
    return f"'{value}" if value[:1] in ("=", "+", "-", "@") else value
