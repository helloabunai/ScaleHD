"""Simulate a sequencing run until we have real data. Named as I remember miseq data runs.
Simulator also provides truth files since we are making the data up.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .seqio import open_text, write_fastq
from .simulate import SimAllele, SimulatedSample, SimulationSpec, simulate
from .structure import AlleleStructure


@dataclass(frozen=True, slots=True)
class RunSample:
    """One sample of a simulated run."""

    name: str
    alleles: tuple[str, ...]
    # Read pairs; None for the run's default.
    pairs: int | None = None
    # Molecules from the last allele, relative to the first.
    abundance: float = 1.0
    # Templates with extra somatic expansion, for alleles of 36 CAG or more.
    somatic_fraction: float = 0.0
    # Only an R1 file, as a single-end run or a lost R2 leaves.
    single_end: bool = False
    # Subfolder of the run it is written in ("" for the run folder itself).
    folder: str = ""


# One of each kind of sample worth looking at, from short normal alleles to ones beyond
# read length, every CCG length seen (from memory), and every atypical structure (from memory).
_FIXED = (
    RunSample("normal-17-21", ("17_1_1_7_2", "21_1_1_7_2")),
    RunSample("homozygous-19", ("19_1_1_7_2",)),
    RunSample("short-9-26", ("9_1_1_7_2", "26_1_1_7_2")),
    RunSample("normal-15-28", ("15_1_1_7_2", "28_1_1_7_2")),
    RunSample("intermediate-17-31", ("17_1_1_7_2", "31_1_1_7_2")),
    RunSample("intermediate-20-35", ("20_1_1_7_2", "35_1_1_7_2")),
    RunSample("reduced-penetrance-18-37", ("18_1_1_7_2", "37_1_1_7_2")),
    RunSample("reduced-penetrance-21-39", ("21_1_1_7_2", "39_1_1_7_2")),
    RunSample("expanded-17-40", ("17_1_1_7_2", "40_1_1_7_2")),
    RunSample("expanded-19-43", ("19_1_1_7_2", "43_1_1_7_2")),
    RunSample("expanded-22-47", ("22_1_1_7_2", "47_1_1_7_2")),
    RunSample("expanded-18-52", ("18_1_1_7_2", "52_1_1_7_2")),
    RunSample("juvenile-20-62", ("20_1_1_7_2", "62_1_1_7_2")),
    RunSample("juvenile-17-70", ("17_1_1_7_2", "70_1_1_7_2")),
    RunSample("juvenile-21-78", ("21_1_1_7_2", "78_1_1_7_2")),
    RunSample("near-read-length-19-82", ("19_1_1_7_2", "82_1_1_7_2")),
    RunSample("beyond-read-length-20-95", ("20_1_1_7_2", "95_1_1_7_2")),
    RunSample("beyond-read-length-18-110", ("18_1_1_7_2", "110_1_1_7_2")),
    RunSample("neighbouring-17-18", ("17_1_1_7_2", "18_1_1_7_2")),
    RunSample("neighbouring-21-22", ("21_1_1_7_2", "22_1_1_7_2")),
    RunSample("ccg-6-12", ("16_1_1_6_2", "24_1_1_12_2")),
    RunSample("ccg-9-10", ("17_1_1_9_2", "43_1_1_10_2")),
    RunSample("ccg-7-10-same-cag", ("17_1_1_7_2", "17_1_1_10_2")),
    RunSample("loss-of-interruption", ("19_1_1_7_2", "42_0_1_7_2")),
    RunSample("caacag-duplication", ("19_2_1_10_2", "40_1_1_7_2")),
    RunSample("ccgcca-deletion", ("19_1_0_7_2", "40_1_1_7_2")),
    RunSample("ccgcca-duplication", ("17_1_1_7_2", "42_1_2_7_2")),
    RunSample("cct-3", ("17_1_1_7_3", "44_1_1_7_2")),
    RunSample("cct-1", ("20_1_1_7_1", "41_1_1_7_2")),
    RunSample("low-depth-17-43", ("17_1_1_7_2", "43_1_1_7_2"), pairs=400),
    RunSample("allele-imbalance-17-43", ("17_1_1_7_2", "43_1_1_7_2"), abundance=0.3),
    RunSample("somatic-expansion-19-44", ("19_1_1_7_2", "44_1_1_7_2"), somatic_fraction=0.2),
    RunSample("r1-only-18-42", ("18_1_1_7_2", "42_1_1_7_2"), single_end=True),
    RunSample("rerun-17-43", ("17_1_1_7_2", "43_1_1_7_2"), folder="reruns"),
    RunSample("rerun-20-45", ("20_1_1_7_2", "45_1_1_7_2"), folder="reruns"),
)

# Normal CAG lengths, most often 17 to 19 (from memory so probably wrong)
_NORMAL = np.arange(9, 27)
_NORMAL_WEIGHTS = np.exp(-0.5 * ((_NORMAL - 18) / 3.0) ** 2)
_NORMAL_WEIGHTS /= _NORMAL_WEIGHTS.sum()
# The second allele (type, share of samples, CAG from, CAG to inclusive).
_SECOND = (
    ("normal", 0.45, 9, 26),
    ("intermediate", 0.10, 27, 35),
    ("reduced", 0.10, 36, 39),
    ("expanded", 0.28, 40, 55),
    ("long", 0.07, 56, 90),
)
_CCG = np.array([7, 10, 9, 6, 12])
_CCG_WEIGHTS = np.array([0.8, 0.12, 0.04, 0.02, 0.02])


def placeholder_samples(random_count: int = 20, seed: int = 1) -> list[RunSample]:
    """The fixed samples above, then ``random_count`` random genotypes (the same per seed).

    Random second alleles are normal, intermediate, reduced penetrance, expanded or long
    in roughly the shares a diagnostic run might see. Expanded alleles mostly have CCG 7,
    and a few have an atypical intervening sequence (from data I have at time of writing).

    Likely will change after feedback/more data/real data
    """
    rng = np.random.default_rng(seed)
    samples = list(_FIXED)
    shares = np.array([share for _, share, _, _ in _SECOND])
    for i in range(random_count):
        first = AlleleStructure(int(rng.choice(_NORMAL, p=_NORMAL_WEIGHTS)))
        first = first.with_counts(ccg=int(rng.choice(_CCG, p=_CCG_WEIGHTS)))
        kind, _, low, high = _SECOND[int(rng.choice(len(_SECOND), p=shares))]
        if kind == "normal":
            second = AlleleStructure(int(rng.choice(_NORMAL, p=_NORMAL_WEIGHTS)))
            second = second.with_counts(ccg=int(rng.choice(_CCG, p=_CCG_WEIGHTS)))
        else:
            second = AlleleStructure(int(rng.integers(low, high + 1)))
            if rng.random() < 0.05:
                second = second.with_counts(ccg=10)
            if second.cag >= 36 and rng.random() < 0.05:
                second = second.with_counts(caacag=0)  # loss of interruption
        if rng.random() < 0.03:
            first = first.with_counts(caacag=2)
        structures = sorted({first, second})
        cags = "-".join(str(x.cag) for x in structures)
        samples.append(RunSample(f"random-{i + 1:02d}-{cags}", tuple(x.label for x in structures)))
    return samples


def write_run(
    folder: str | Path,
    samples: Sequence[RunSample],
    *,
    pairs: int = 20_000,
    seed: int = 1,
    workers: int | None = None,
) -> list[Path]:
    """Write ``samples`` as a MiSeq run folder. Returns every file written.

    Refuses a folder that already has files in it, so a run is never mixed with another.
    ``workers`` samples are simulated at once (default is every core).
    """
    folder = Path(folder)
    if folder.exists() and any(folder.iterdir()):
        raise FileExistsError(f"{folder} already has files in it")
    tasks = [
        _Task(folder / s.folder, f"{s.name}_S{n}_L001", s, s.pairs or pairs, seed * 1000 + n)
        for n, s in enumerate(samples, start=1)
    ]
    if workers == 1:
        written = [path for task in tasks for path in _write_sample(task)]
    else:
        with ProcessPoolExecutor(workers) as pool:
            written = [path for paths in pool.map(_write_sample, tasks) for path in paths]

    # Reads that matched no sample's index: junk, not the repeat.
    junk = simulate(
        SimulationSpec(
            (SimAllele(AlleleStructure(17)),), pairs=max(pairs // 20, 10), off_target=1.0, seed=seed
        )
    )
    written += _write_reads(folder, "Undetermined_S0_L001", junk, single_end=False)
    # An R2 whose R1 is missing, as a copy gone wrong leaves.
    orphan = simulate(
        SimulationSpec(
            (SimAllele(AlleleStructure(19)), SimAllele(AlleleStructure(44))),
            pairs=max(pairs // 10, 10),
            seed=seed * 1000 + len(samples) + 1,
        )
    )
    path = folder / f"orphan_S{len(samples) + 1}_L001_R2_001.fastq.gz"
    with open_text(path, "wt") as handle:
        write_fastq(handle, orphan.r2)
    written.append(path)

    readme = folder / "README.txt"
    readme.write_text(
        "Generated placeholder data, not real sequencing.\n\n"
        f"Written by `scalehd simulate-run` (seed {seed}, {pairs:,} read pairs per sample\n"
        "unless a sample's truth file says otherwise) to build and test ScaleHD's job\n"
        "pipeline before real data is available.\n"
    )
    written.append(readme)
    return written


@dataclass(frozen=True, slots=True)
class _Task:
    folder: Path
    prefix: str  # file name up to _R1_001.fastq.gz
    sample: RunSample
    pairs: int
    seed: int


def _write_sample(task: _Task) -> list[Path]:
    s = task.sample
    structures = [AlleleStructure.from_label(label) for label in s.alleles]
    alleles = [
        SimAllele(x, somatic_fraction=s.somatic_fraction if x.cag >= 36 else 0.0)
        for x in structures
    ]
    alleles[-1] = SimAllele(alleles[-1].structure, s.abundance, alleles[-1].somatic_fraction)
    simulated = simulate(SimulationSpec(tuple(alleles), pairs=task.pairs, seed=task.seed))
    written = _write_reads(task.folder, task.prefix, simulated, single_end=s.single_end)
    truth = task.folder / f"{s.name}.truth.json"
    truth.write_text(json.dumps(simulated.truth(s.name), indent=2) + "\n")
    return [*written, truth]


def _write_reads(
    folder: Path, prefix: str, simulated: SimulatedSample, *, single_end: bool
) -> list[Path]:
    folder.mkdir(parents=True, exist_ok=True)
    written = []
    for read, records in ((1, simulated.r1), (2, simulated.r2)):
        if read == 2 and single_end:
            continue
        path = folder / f"{prefix}_R{read}_001.fastq.gz"
        with open_text(path, "wt") as handle:
            write_fastq(handle, records)
        written.append(path)
    return written
