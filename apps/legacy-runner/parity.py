"""Automated process to confirm the legacy wrapper gets the same exact output as ScaleHD 1.x does
on its own, for simulated samples. First go is all as a batch via CLI ScaleHD. Then, run one
at a time through the wrapper (as would happen when called by the web frontend).
All of it runs two times, once with demultiplexing and once without. Also run
the lab's ContinuousTesting sample file through to verify (same test file as from 1.3 legacy).

Before any of that, contract.py checks 1.x still has what the new code relies on.

    python apps/legacy-runner/parity.py [--work DIR]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any
from xml.etree import ElementTree

import contract
import scalehd_legacy_run as runner
from scalehd.simulate import SimAllele, SimulationSpec, simulate
from scalehd.structure import AlleleStructure

LEGACY_ROOT = Path("/opt/legacy")
RUNNER = LEGACY_ROOT / "bin" / "scalehd-legacy-run"
SCALEHD_1X = LEGACY_ROOT / "python" / "bin" / "ScaleHD"
LAB_TEST = contract.UPSTREAM / "ContinuousTesting"

# label: (alleles, read pairs)
SAMPLES = {
    "normal": (["17_1_1_7_2", "21_1_1_7_2"], 4000),
    "expanded": (["19_1_1_7_2", "44_1_1_7_2"], 4000),
    "homozygous": (["18_1_1_7_2"], 4000),
    "ccg": (["17_1_1_7_2", "20_1_1_10_2"], 4000),
    "atypical": (["17_1_1_7_2", "22_2_1_7_2"], 4000),
    "few-reads": (["17_1_1_7_2", "21_1_1_7_2"], 150),
}

# 1.x's documented configuration (docs/ProvidingInput.rst)
CONFIG = {
    "instance_flags": {
        "demultiplex": "False",
        "quality_control": "True",
        "sequence_alignment": "True",
        "atypical_realignment": "True",
        "genotype_prediction": "True",
        "snp_calling": "True",
    },
    "demultiplex_flags": {
        "forward_adapter": "GCGACCCTGG",
        "forward_position": "5P",
        "reverse_adapter": "GCAGCGGCTG",
        "reverse_position": "5P",
        "error_rate": "0",
        "min_overlap": "10",
        "min_length": "",
        "max_length": "",
    },
    "trim_flags": {
        "trim_type": "Adapter",
        "quality_threshold": "5",
        "adapter_flag": "-a",
        "forward_adapter": "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        "reverse_adapter": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        "error_tolerance": "0.39",
    },
    "alignment_flags": {
        "min_seed_length": "19",
        "band_width": "100",
        "seed_length_extension": "1.5",
        "skip_seed_with_occurrence": "500",
        "chain_drop": "0.50",
        "seeded_chain_drop": "0",
        "seq_match_score": "1",
        "mismatch_penalty": "4",
        "indel_penalty": "6,6",
        "gap_extend_penalty": "6,6",
        "prime_clipping_penalty": "5,5",
        "unpaired_pairing_penalty": "17",
    },
    "prediction_flags": {"snp_observation_threshold": "2", "quality_cutoff": "0"},
}
DEMULTIPLEXING = {
    **CONFIG,
    "instance_flags": {**CONFIG["instance_flags"], "demultiplex": "True"},
}
CONFIGS = {
    "documented": ("1.x's documented settings", CONFIG),
    "demultiplexing": ("the same, demultiplexing first", DEMULTIPLEXING),
}


def simulate_samples(folder: Path) -> dict[str, tuple[Path, Path]]:
    reads = {}
    for label, (alleles, pairs) in SAMPLES.items():
        spec = SimulationSpec(
            tuple(SimAllele(AlleleStructure.from_label(allele)) for allele in alleles),
            pairs=pairs,
        )
        r1, r2, _ = simulate(spec).write(folder, label)
        reads[label] = (r1, r2)
    return reads


def run_batch(
    reads: dict[str, tuple[Path, Path]], config: dict[str, dict[str, str]], work: Path
) -> Path:
    """Run ScaleHD 1.x as if on CLI"""
    data = work / "data"
    data.mkdir(parents=True)
    for r1, r2 in reads.values():
        shutil.copy(r1, data)
        shutil.copy(r2, data)
    config_file = work / "config.xml"
    references = LEGACY_ROOT / "references"
    config_file.write_text(
        runner.config_xml(
            config,
            str(data),
            str(references / runner.FORWARD_REFERENCE),
            str(references / runner.REVERSE_REFERENCE),
        )
    )
    paths = [str(LEGACY_ROOT / folder) for folder in ("bin", "python/bin", "tools/bin")]
    env = {**os.environ, "PATH": ":".join([*paths, os.environ["PATH"]])}
    with open(work / "batch.log", "w") as log:
        subprocess.run(
            [str(SCALEHD_1X), "-v", "-c", str(config_file), "-o", str(work), "-j", "batch"],
            env=env,
            stdin=subprocess.DEVNULL,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    return work / "batch"


def run_sample(
    label: str, r1: Path, r2: Path, config: dict[str, dict[str, str]], threads: int, work: Path
) -> dict[str, Any]:
    request = work / f"{label}.request.json"
    request.write_text(
        json.dumps(
            {"label": label, "r1": str(r1), "r2": str(r2), "threads": threads, "config": config}
        )
    )
    out = work / label
    with open(work / f"{label}.log", "w") as log:
        done = subprocess.run(
            [str(RUNNER), str(request), str(out)],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=log,
            text=True,
        )
    if done.returncode not in (0, 1):
        raise SystemExit(f"the runner failed on {label} ({threads} threads): see {work}")
    result: dict[str, Any] = json.loads((out / "result.json").read_text())
    stages = [json.loads(line)["stage"] for line in done.stdout.splitlines()]
    if stages != result["stages"]:
        raise SystemExit(f"{label}: stdout said {stages}, result.json {result['stages']}")
    return _without(result, str(out))


def _without(value: Any, out: str) -> Any:
    """result.json without what differs between any two runs i.e. its folder, timing, threads."""
    if isinstance(value, dict):
        return {k: _without(v, out) for k, v in value.items() if k not in ("seconds", "threads")}
    if isinstance(value, list):
        return [_without(item, out) for item in value]
    if isinstance(value, str):
        return value.replace(out, "<out>")
    return value


def rows(text: str | None) -> dict[str, list[str]]:
    """A 1.x report's rows by sample label (its first column)"""
    found: dict[str, list[str]] = {}
    for line in (text or "").splitlines()[1:]:
        found.setdefault(line.split(",", 1)[0], []).append(line)
    return found


def differences(a: Any, b: Any, where: str = "") -> list[str]:
    if isinstance(a, dict) and isinstance(b, dict):
        return [
            difference
            for key in sorted(set(a) | set(b))
            for difference in differences(a.get(key), b.get(key), f"{where}.{key}")
        ]
    return [] if a == b else [f"{where}: {a!r} != {b!r}"[:300]]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--work", type=Path, help="keep everything here (default: a temp folder)")
    args = parser.parse_args()
    work = args.work or Path(tempfile.mkdtemp(prefix="legacy-parity-"))
    work.mkdir(parents=True, exist_ok=True)

    print(f"What the runner and server rely on in 1.x ({contract.UPSTREAM})", flush=True)
    failures = [f"contract, {failure}" for failure in contract.check()]

    print(f"\nSimulating {len(SAMPLES)} samples in {work}", flush=True)
    reads = simulate_samples(work / "reads")
    for folder, (title, config) in CONFIGS.items():
        print(f"\nWith {title}", flush=True)
        failures += [f"{folder}, {failure}" for failure in check(reads, config, work / folder)]

    print("\nThe Monckton lab's own test sample", flush=True)
    failures += [f"lab test, {failure}" for failure in check_lab_test(work / "lab-test")]

    print()
    for failure in failures:
        print(failure)
    print("Parity! yay!" if not failures else f"{len(failures)} differences.")
    return 1 if failures else 0


def check_lab_test(work: Path) -> list[str]:
    """The lab's sample, configuration and references, run as their test runs 1.x (with -b):
    the reports must match the lab's expected SHA-256s."""
    config_file = LAB_TEST / "ConfigurationFile_MiSeq_runX.xml"
    root = ElementTree.fromstring(config_file.read_text())
    (r1,) = (LAB_TEST / "FastqFiles").glob("*_R1.fastq.gz")
    request = {
        # 1.x's name for the sample, from its files.
        "label": r1.name.split(".")[0].removesuffix("_R1"),
        "r1": str(r1),
        "r2": str(r1.with_name(r1.name.replace("_R1", "_R2"))),
        "forward_reference": str(
            LAB_TEST / "References" / Path(root.attrib["forward_reference"]).name
        ),
        "reverse_reference": str(
            LAB_TEST / "References" / Path(root.attrib["reverse_reference"]).name
        ),
        "threads": max(3, os.cpu_count() or 1),
        "broadscope": True,
        "config": {element.tag: dict(element.attrib) for element in root},
    }
    work.mkdir(parents=True)
    (work / "request.json").write_text(json.dumps(request))
    with open(work / "runner.log", "w") as log:
        subprocess.run(
            [str(RUNNER), str(work / "request.json"), str(work / "out")],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=log,
        )
    expected = dict(
        line.split(",")
        for line in (LAB_TEST / "scalehd_test_expected_checksums.csv").read_text().split()
    )
    failures = []
    for report, digest in expected.items():
        path = work / "out" / runner.JOB_NAME / report
        found = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        print(f"{report:<28} {'ok' if found == digest else 'DIFFERS'}")
        if found != digest:
            failures.append(f"{report} isn't what the lab expects")
    return failures


def check(
    reads: dict[str, tuple[Path, Path]], config: dict[str, dict[str, str]], work: Path
) -> list[str]:
    """Runs the samples as a batch and through the runner, prints how each went, and returns
    what differed."""
    many = max(3, os.cpu_count() or 1)
    print("ScaleHD 1.x on its own, as one batch", flush=True)
    batch = run_batch(reads, config, work / "batch")
    batch_report = rows((batch / "InstanceReport.csv").read_text())
    batch_distributions = rows((batch / "AlignedDistributions.csv").read_text())

    print("The runner w/ all samples at once, 1 thread each", flush=True)
    (work / "one").mkdir()
    with ThreadPoolExecutor(len(SAMPLES)) as pool:
        running = {
            label: pool.submit(run_sample, label, *pair, config, 1, work / "one")
            for label, pair in reads.items()
        }
        one = {label: future.result() for label, future in running.items()}
    print(f"The runner w/ one sample at a time, {many} threads each", flush=True)
    (work / "many").mkdir()
    several = {
        label: run_sample(label, *pair, config, many, work / "many")
        for label, pair in reads.items()
    }

    failures = []
    print()
    print(f"{'sample':<12} {'simulated':<24} {'1.x genotype':<24} {'status':<16} parity")
    for label, (alleles, _) in SAMPLES.items():
        result = one[label]
        sample = result["sample"] or {}
        # 1.x's name for the sample, DMPX_<label> once demultiplexed.
        name = sample.get("label", label)
        problems = []
        if rows(result["instance_report"]) != {name: batch_report.get(name)}:
            problems.append("InstanceReport row differs from the batch")
        if rows(result["aligned_distributions"]).get(name) != batch_distributions.get(name):
            problems.append("AlignedDistributions rows differ from the batch")
        changed = differences(result, several[label])
        if changed:
            problems.append(f"{len(changed)} values differ with {many} threads: {changed[:3]}")
        genotype = "/".join(
            (sample.get(allele) or {}).get("allelegenotype") or "-"
            for allele in ("primaryallele", "secondaryallele")
        )
        status = result["status"] + (
            f" ({result['failed_stage']})" if result["failed_stage"] else ""
        )
        simulated = "/".join(alleles if len(alleles) == 2 else alleles * 2)
        parity = "DIFFERS" if problems else "ok"
        print(f"{label:<12} {simulated:<24} {genotype:<24} {status:<16} {parity}")
        failures += [f"{label}: {problem}" for problem in problems]
    return failures


if __name__ == "__main__":
    sys.exit(main())
