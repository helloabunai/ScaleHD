"""ScaleHD legacy genotyping. Version 1.x ran through a wrapper, and reading output.

The wrapper (apps/legacy-runner/scalehd_legacy_run.py, in the Docker image) runs 1.x on one
sample, unchanged, and saves everything 1.x knew about it in result.json. In the sample's
folder:

    legacy-request.json   what the runner was asked to do
    legacy.log            1.x's log, and the output of the programs it ran
    legacy/               output - generated config.xml, run/ (raw results), result.json
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

from scalehd.structure import FIELDS, AlleleStructure

from .progress import ProgressReporter, SampleStage
from .schemas import LegacySettings

FOLDER = "legacy"
REQUEST = "legacy-request.json"
LOG = "legacy.log"
RESULT = "result.json"
CALL_SCHEMA = "scalehd.legacy-call/1"
# A sample name 1.x can use, as the runner checks it.
_LABEL = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]*")
# 1.x's reports for the whole job i.e. batch only
REPORTS = {
    "InstanceReport": "instance_report",
    "AlignedDistributions": "aligned_distributions",
}

# 1.x's name for the stage a sample failed in (its report's Exception Raised)
FAILED_STAGES = {
    "SeqQC": "quality control",
    "SeqALN": "alignment",
    "DSP": "atypical scanning",
    "SeqRE-ALN": "realignment",
    "Genotype": "genotyping",
    "SNPCalling": "SNP calling",
    "Report/Graph": "writing its report",
}
# 1.x's sample flags: the InstanceReport column, and the value it comes from in result.json.
FLAGS = (
    ("Homozygous Haplotype", "homozygoushaplotype"),
    ("Neighbouring Peaks", "neighbouringpeaks"),
    ("Diminished Peaks", "diminishedpeaks"),
    ("Novel Atypical", "novel_atypical_structure"),
    ("Alignment Warning", "alignmentwarning"),
    ("Atypical Alignment Warning", "atypical_alignmentwarning"),
    ("CCG Rewritten", "atypical_ccgrewrite"),
    ("CCG Zygosity Rewritten", "atypical_zygrewrite"),
    ("CCG Uncertainty", "ccguncertainty"),
    ("CCT Uncertainty", "cctuncertainty"),
    ("SVM Failure", "svm_failure"),
    ("Differential Confusion", "differential_confusion"),
    ("Missed Expansion", "missed_expansion"),
    ("Heuristic Filtering Success", "heuristicfilter"),
    ("Peak Inspection Warning", "peakinspection_warning"),
    ("Low Distribution Reads", "distribution_readcount_warning"),
    ("Low Peak Reads", "fatalreadallele"),
)
# The one flag that says something went right (lol) (1.x's HTML report called it FILTERPASS).
PASS_FLAG = "Heuristic Filtering Success"
# 1.x's forward distribution is aligned reads per reference, CAG 1-200 for each CCG in turn.
CAG_LENGTHS = 200
# How long each 1.x stage takes, as parts of a sample's run
# Timed on simulated and real data for rough progress bar feedback from the wrapper to UI
STAGE_PARTS = {
    SampleStage.PREPARING: (1.1, 0.0),
    SampleStage.INDEXING: (0.2, 0.0),
    SampleStage.DEMULTIPLEXING: (0.3, 0.0),
    SampleStage.QUALITY_CONTROL: (2.7, 0.0),
    SampleStage.ALIGNMENT: (31.9, 31.7),
    SampleStage.ATYPICAL_SCANNING: (9.6, 7.2),
    SampleStage.GENOTYPING: (4.3, 0.0),
    SampleStage.SNP_CALLING: (1.9, 0.0),
}


class LegacyError(Exception):
    """1.x couldn't genotype the sample, or the runner couldn't run it."""


def label(name: str, folder: Path) -> str:
    return name if _LABEL.fullmatch(name) else folder.name


def run_sample(
    name: str,
    folder: Path,
    r1: Path,
    r2: Path | None,
    settings: LegacySettings,
    threads: int,
    runner: Path,
    progress: ProgressReporter,
) -> dict[str, Any]:
    """Run 1.x on one sample. Returns result.json, or raises LegacyError"""
    if r2 is None:
        raise LegacyError("Legacy ScaleHD needs both read files (R1 and R2)")
    out_dir = folder / FOLDER
    # Left by an earlier attempt: the runner needs an empty folder.
    shutil.rmtree(out_dir, ignore_errors=True)
    folder.mkdir(parents=True, exist_ok=True)
    request = {
        "label": label(name, folder),
        "r1": str(r1),
        "r2": str(r2),
        "threads": threads,
        "enshrine": settings.enshrine,
        "broadscope": settings.broadscope,
        "config": settings.config(),
    }
    (folder / REQUEST).write_text(json.dumps(request, indent=2) + "\n")
    done = stage_progress(settings, threads)
    with open(folder / LOG, "w") as log:
        process = subprocess.Popen(
            [str(runner), str(folder / REQUEST), str(out_dir)],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=log,
            text=True,
        )
        assert process.stdout is not None
        for line in process.stdout:
            stage = _stage(line)
            if stage is not None:
                progress.stage(stage, done.get(stage))
        code = process.wait()
    result = read_result(folder)
    if result is None:
        raise LegacyError(f"the legacy runner stopped (exit {code}): {_last_line(folder / LOG)}")
    if result.get("status") != "complete":
        raise LegacyError(failure(result))
    return result


def stage_progress(settings: LegacySettings, cores: int = 1) -> dict[SampleStage, float]:
    """How much of a sample's run is done as each stage starts, for the stages these
    settings run, on this many cores."""
    flags = settings.instance_flags
    left_out = {
        SampleStage.DEMULTIPLEXING: not flags.demultiplex,
        SampleStage.QUALITY_CONTROL: not flags.quality_control,
        SampleStage.SNP_CALLING: not flags.snp_calling,
    }
    stages = [stage for stage in STAGE_PARTS if not left_out.get(stage, False)]
    seconds = {stage: stage_seconds(stage, cores) for stage in stages}
    total = sum(seconds.values())
    done: dict[SampleStage, float] = {}
    so_far = 0.0
    for stage in stages:
        done[stage] = so_far / total
        so_far += seconds[stage]
    # Between atypical scanning and genotyping.
    done[SampleStage.REALIGNMENT] = done[SampleStage.GENOTYPING]
    return done


def stage_seconds(stage: SampleStage, cores: int) -> float:
    """Roughly how long a stage takes on this many cores."""
    one_core, shared = STAGE_PARTS[stage]
    shared_by = cores
    if stage is SampleStage.ATYPICAL_SCANNING:
        shared_by = 3 if cores >= 3 else 1
    return one_core - shared + shared / shared_by


def _stage(line: str) -> SampleStage | None:
    try:
        return SampleStage(json.loads(line)["stage"])
    except (ValueError, KeyError, TypeError):
        return None


def _last_line(path: Path) -> str:
    try:
        lines = path.read_text(errors="replace").strip().splitlines()
    except OSError:
        return "no log"
    return lines[-1] if lines else "no log"


def read_result(folder: Path) -> dict[str, Any] | None:
    try:
        result: dict[str, Any] = json.loads((folder / FOLDER / RESULT).read_text())
    except (OSError, ValueError):
        return None
    return result


def failure(result: dict[str, Any]) -> str:
    stage = result.get("failed_stage")
    where = f" in {FAILED_STAGES.get(stage, stage)}" if stage else ""
    return f"Legacy ScaleHD failed{where}: {result.get('error') or 'it gave no reason'}"


def call(result: dict[str, Any]) -> dict[str, Any]:
    sample = result["sample"]
    alleles = [_allele(sample[key]) for key in ("primaryallele", "secondaryallele")]
    values = {column: bool(sample.get(key)) for column, key in FLAGS}
    return {
        "schema": CALL_SCHEMA,
        # Primary then secondary, as 1.x orders them.
        "genotype": "/".join(allele["structure"] for allele in alleles),
        "alleles": alleles,
        # The flags 1.x raised, and every flag's value.
        "flags": [column for column, raised in values.items() if raised and column != PASS_FLAG],
        "flag_values": values,
        "versions": result.get("versions"),
    }


def _allele(values: dict[str, Any]) -> dict[str, Any]:
    # 1.x's final label for the allele (its report's GTYPE column).
    structure = values["reflabel"]
    counts: dict[str, int | None]
    try:
        counts = dict(zip(FIELDS, AlleleStructure.from_label(structure).counts, strict=True))
    except ValueError:
        counts = dict.fromkeys(FIELDS)
    return {
        "structure": structure,
        "status": values["allelestatus"],
        **counts,
        "intervening": values["intervening"],
        "confidence": _number(values["alleleconfidence"]),
        "confidence_interval": values["alleleconfinterval"],
        "backward_slippage": _number(values["backwardsslippage"]),
        "somatic_mosaicism": _number(values["somaticmosaicism"]),
        "variant_call": values["variantcall"],
        "variant_score": _number(values["variantscore"]),
        "map_fw": _number(values["fwalncount"]),
        "map_percent_fw": _number(values["fwalnpcnt"]),
        "purged_fw": _number(values["fwalnrmvd"]),
        "map_rv": _number(values["rvalncount"]),
        "map_percent_rv": _number(values["rvalnpcnt"]),
        "purged_rv": _number(values["rvalnrmvd"]),
    }


def _number(value: Any) -> Any:
    """A number 1.x kept as text (it scrapes some from samtools' output) as a number."""
    if isinstance(value, str):
        for kind in (int, float):
            try:
                return kind(value)
            except ValueError:
                pass
    return value


def confidence(legacy_call: dict[str, Any]) -> float | None:
    """The sample's confidence, for job lists 0-100."""
    values = [allele["confidence"] for allele in legacy_call["alleles"]]
    numbers = [value for value in values if isinstance(value, int | float)]
    return float(min(numbers)) if numbers else None


def matches(legacy_call: dict[str, Any], truth: str) -> bool | None:
    """Whether 1.x's genotype is the true one, ignoring allele order."""
    try:
        called = sorted(AlleleStructure.from_label(a["structure"]) for a in legacy_call["alleles"])
    except ValueError:
        return False
    return called == sorted(AlleleStructure.from_label(label) for label in truth.split("/"))


def cag_reads(allele: dict[str, Any]) -> list[int]:
    """Aligned reads at each CAG length (1-200) at the allele's CCG, from the forward reads"""
    start = (allele["ccg"] - 1) * CAG_LENGTHS
    reads: list[int] = allele["fwarray"][start : start + CAG_LENGTHS]
    return reads


def graphs(result: dict[str, Any]) -> Path | None:
    """1.x's PDF of the sample's graphs."""
    sample = result.get("sample") or {}
    if not sample.get("predictpath") or not sample.get("label"):
        return None
    return Path(sample["predictpath"]) / f"{sample['label']}.pdf"


def fastqc(result: dict[str, Any]) -> Path | None:
    """FastQC's own HTML report on the forward reads."""
    reports = (result.get("sample") or {}).get("fqcreport") or []
    return Path(reports[0]).with_name("fastqc_report.html") if reports else None


def job_report(folders: list[Path], report: str) -> str | None:
    """One of 1.x's reports for the whole job."""
    key = REPORTS[report]
    header: str | None = None
    rows: list[str] = []
    for folder in folders:
        result = read_result(folder)
        text = result.get(key) if result else None
        if not text:
            continue
        first, *lines = text.splitlines(keepends=True)
        header = header or first
        rows += lines
    return None if header is None else header + "".join(rows)
