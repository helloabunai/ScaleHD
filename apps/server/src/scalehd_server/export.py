"""A job's results as one CSV table: a row per sample, its two alleles side by side.

Like ScaleHD 1.x's InstanceReport.csv, for opening in a spreadsheet. Built from the calls
stored with each sample, so samples not yet called have their results left blank.
"""

from __future__ import annotations

import csv
import io
from typing import Any

from .models import Job, Sample
from .workspace import folder_name

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


def results_filename(job: Job) -> str:
    label = "demo" if job.demo else folder_name(job.name, "job")
    return f"{job.id}-{label}-results.csv"


def results_csv(job: Job) -> str:
    out = io.StringIO()
    writer = csv.writer(out, lineterminator="\n")
    writer.writerow(HEADER)
    for sample in job.samples:
        writer.writerow(_row(sample))
    return out.getvalue()


def _row(sample: Sample) -> list[str]:
    call = sample.call or {}
    row = [
        _text(sample.name),
        str(sample.status),
        call.get("genotype", ""),
        _number(call.get("quality"), 1),
        _number(call.get("posterior"), 4),
        "; ".join(call.get("flags", [])),
        _number(call.get("molecules"), 0),
    ]
    alleles: list[dict[str, Any]] = call.get("alleles", [])
    for n in (0, 1):
        allele = alleles[n] if n < len(alleles) else {}
        row.append(allele.get("structure", ""))
        for _, key in _ALLELE:
            row.append(_allele_value(allele, key))
    row += [
        sample.truth or "",
        {True: "yes", False: "no", None: ""}[sample.matches_truth],
        _text(sample.error or ""),
        _text(sample.r1 or ""),
        _text(sample.r2 or ""),
    ]
    return row


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


def _number(value: float | None, places: int) -> str:
    if value is None:
        return ""
    return f"{value:.{places}f}"


def _text(value: str) -> str:
    """Free text, kept from being read as a formula by a spreadsheet."""
    return f"'{value}" if value[:1] in ("=", "+", "-", "@") else value
