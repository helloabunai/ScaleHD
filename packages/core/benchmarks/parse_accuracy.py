"""Per-molecule parsing accuracy on simulated samples.

A joined molecule is wrong if an exact field differs from the simulated truth or a
lower bound exceeds it. Run with ``uv run python packages/core/benchmarks/parse_accuracy.py``.

Really, truly not exhaustive as a benchmark. Consider WIP.
"""

from __future__ import annotations

import argparse
import time
from collections import Counter

from scalehd.pairs import join_mates
from scalehd.parse import RepeatParser
from scalehd.seqio import reverse_complement
from scalehd.simulate import SimAllele, SimulationSpec, simulate
from scalehd.structure import AlleleStructure, FieldStatus, Observation

SCENARIOS = [
    ("17_1_1_7_2", "43_1_1_7_2"),
    ("42_0_1_7_2", "19_2_1_10_2"),
    ("21_1_1_7_2",),
    ("20_1_1_7_2", "66_1_1_7_2"),
    ("20_1_1_7_2", "75_1_1_7_2"),
    ("20_1_1_7_2", "95_1_1_7_2"),
    ("18_1_1_7_2", "120_1_1_7_2"),
]


def consistent(observation: Observation, truth: AlleleStructure) -> bool:
    for seen, status, actual in zip(
        observation.counts, observation.status, truth.counts, strict=True
    ):
        if status is FieldStatus.EXACT and seen != actual:
            return False
        if status is FieldStatus.LOWER_BOUND and seen > actual:
            return False
    return True


def run(labels: tuple[str, ...], pairs: int, seed: int) -> tuple[Counter[str], float]:
    alleles = tuple(
        SimAllele(s, somatic_fraction=0.1 if s.cag > 35 else 0.0)
        for s in map(AlleleStructure.from_label, labels)
    )
    sample = simulate(SimulationSpec(alleles, pairs=pairs, seed=seed))
    parser = RepeatParser()
    tally: Counter[str] = Counter()
    started = time.perf_counter()
    for r1, r2 in zip(sample.r1, sample.r2, strict=True):
        a = parser.parse(r1.sequence)
        b = parser.parse(reverse_complement(r2.sequence))
        joined = join_mates(
            a.observation if a.usable else None, b.observation if b.usable else None
        )
        observation = joined.observation
        if observation is None or observation.is_empty:
            tally["dropped" if joined.is_discordant else "unusable"] += 1
            continue
        truth = AlleleStructure.from_label(r1.name.rsplit(":", 1)[-1])
        tally["complete" if observation.is_complete else "partial"] += 1
        tally["wrong"] += not consistent(observation, truth)
    return tally, time.perf_counter() - started


def main() -> None:
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("--pairs", type=int, default=10_000)
    args.add_argument("--seed", type=int, default=11)
    options = args.parse_args()
    print(f"{'genotype':<26}{'complete':>9}{'partial':>9}{'dropped':>9}{'wrong':>7}{'secs':>7}")
    for labels in SCENARIOS:
        tally, seconds = run(labels, options.pairs, options.seed)
        print(
            f"{'/'.join(labels):<26}{tally['complete']:>9}{tally['partial']:>9}"
            f"{tally['dropped']:>9}{tally['wrong']:>7}{seconds:>7.1f}"
        )


if __name__ == "__main__":
    main()
