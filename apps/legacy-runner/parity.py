"""Automated process to confirm the legacy wrapper gets the same exact output as ScaleHD 1.x does
on its own, for simulated samples. First go is all as a batch via CLI ScaleHD. Then, run one
at a time through the wrapper (as would happen when called by the web frontend).

    python apps/legacy-runner/parity.py [--work DIR]
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import scalehd_legacy_run as runner
from scalehd.simulate import SimAllele, SimulationSpec, simulate
from scalehd.structure import AlleleStructure

LEGACY_ROOT = Path("/opt/legacy")
RUNNER = LEGACY_ROOT / "bin" / "scalehd-legacy-run"
SCALEHD_1X = LEGACY_ROOT / "python" / "bin" / "ScaleHD"

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


def run_batch(reads: dict[str, tuple[Path, Path]], work: Path) -> Path:
    """Run ScaleHD 1.x as if on CLI"""
    data = work / "data"
    data.mkdir(parents=True)
    for r1, r2 in reads.values():
        shutil.copy(r1, data)
        shutil.copy(r2, data)
    config = work / "config.xml"
    references = LEGACY_ROOT / "references"
    config.write_text(
        runner.config_xml(
            CONFIG,
            str(data),
            str(references / runner.FORWARD_REFERENCE),
            str(references / runner.REVERSE_REFERENCE),
        )
    )
    paths = [str(LEGACY_ROOT / folder) for folder in ("bin", "python/bin", "tools/bin")]
    env = {**os.environ, "PATH": ":".join([*paths, os.environ["PATH"]])}
    with open(work / "batch.log", "w") as log:
        subprocess.run(
            [str(SCALEHD_1X), "-v", "-c", str(config), "-o", str(work), "-j", "batch"],
            env=env,
            stdin=subprocess.DEVNULL,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    return work / "batch"


def run_sample(label: str, r1: Path, r2: Path, threads: int, work: Path) -> dict[str, Any]:
    request = work / f"{label}.request.json"
    request.write_text(
        json.dumps(
            {"label": label, "r1": str(r1), "r2": str(r2), "threads": threads, "config": CONFIG}
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
    many = max(3, os.cpu_count() or 1)

    print(f"Simulating {len(SAMPLES)} samples in {work}", flush=True)
    reads = simulate_samples(work / "reads")
    print("ScaleHD 1.x on its own, as one batch", flush=True)
    batch = run_batch(reads, work / "batch")
    batch_report = rows((batch / "InstanceReport.csv").read_text())
    batch_distributions = rows((batch / "AlignedDistributions.csv").read_text())

    print("The runner w/ all samples at once, 1 thread each", flush=True)
    (work / "one").mkdir()
    with ThreadPoolExecutor(len(SAMPLES)) as pool:
        running = {
            label: pool.submit(run_sample, label, *pair, 1, work / "one")
            for label, pair in reads.items()
        }
        one = {label: future.result() for label, future in running.items()}
    print(f"The runner w/ one sample at a time, {many} threads each", flush=True)
    (work / "many").mkdir()
    several = {
        label: run_sample(label, *pair, many, work / "many") for label, pair in reads.items()
    }

    failures = []
    print()
    print(f"{'sample':<12} {'simulated':<24} {'1.x genotype':<24} {'status':<16} parity")
    for label, (alleles, _) in SAMPLES.items():
        result = one[label]
        problems = []
        if rows(result["instance_report"]) != {label: batch_report.get(label)}:
            problems.append("InstanceReport row differs from the batch")
        if rows(result["aligned_distributions"]).get(label) != batch_distributions.get(label):
            problems.append("AlignedDistributions rows differ from the batch")
        changed = differences(result, several[label])
        if changed:
            problems.append(f"{len(changed)} values differ with {many} threads: {changed[:3]}")
        sample = result["sample"] or {}
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

    print()
    for failure in failures:
        print(failure)
    print("Parity! yay!" if not failures else f"{len(failures)} differences.")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
