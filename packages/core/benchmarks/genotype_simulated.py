"""Genotype-call accuracy on simulated samples.

Each scenario is simulated with several seeds and called from FASTQ through parsing,
read pair joining and the caller. A call is right when its genotype label matches the truth,
for an allele beyond read length, when the label is a lower bound at or below the
true CAG. Run with ``uv run python packages/core/benchmarks/genotype_simulated.py``.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

from scalehd.counts import count_reads
from scalehd.genotype import GenotypeCall, call_genotype
from scalehd.simulate import SimAllele, SimulationSpec, call_matches, simulate, true_genotype
from scalehd.structure import AlleleStructure


@dataclass(frozen=True)
class Scenario:
    name: str
    alleles: tuple[str, ...]
    pairs: int = 5000
    somatic: float = 0.0  # extra somatic expansion for alleles of 36 CAG or more
    abundance: float = 1.0  # of the last allele


SCENARIOS = (
    Scenario("normal heterozygote", ("17_1_1_7_2", "21_1_1_7_2")),
    Scenario("expanded", ("17_1_1_7_2", "43_1_1_7_2")),
    Scenario("expanded, extra somatic", ("19_1_1_7_2", "44_1_1_7_2"), somatic=0.2),
    Scenario("homozygous", ("21_1_1_7_2",)),
    Scenario("neighbouring +1", ("17_1_1_7_2", "18_1_1_7_2")),
    Scenario("neighbouring -1", ("16_1_1_7_2", "17_1_1_7_2")),
    Scenario("same CAG, CCG 7/10", ("17_1_1_7_2", "17_1_1_10_2")),
    Scenario("loss of interruption", ("42_0_1_7_2", "19_1_1_7_2")),
    Scenario("CAACAG duplication", ("19_2_1_10_2", "40_1_1_7_2")),
    Scenario("long, 55", ("20_1_1_7_2", "55_1_1_7_2")),
    Scenario("long, 60", ("20_1_1_7_2", "60_1_1_7_2")),
    Scenario("long", ("20_1_1_7_2", "66_1_1_7_2")),
    Scenario("long, 70", ("20_1_1_7_2", "70_1_1_7_2")),
    Scenario("very long", ("20_1_1_7_2", "75_1_1_7_2")),
    # Around read length. some molecules are read in full, the rest only to a lower call boundary.
    Scenario("at read length, 80", ("20_1_1_7_2", "80_1_1_7_2")),
    Scenario("at read length, 84", ("20_1_1_7_2", "84_1_1_7_2")),
    Scenario("beyond read length", ("20_1_1_7_2", "95_1_1_7_2")),
    Scenario("low depth", ("17_1_1_7_2", "43_1_1_7_2"), pairs=300),
    Scenario("allele imbalance", ("17_1_1_7_2", "43_1_1_7_2"), abundance=0.3),
)


def run(task: tuple[Scenario, int]) -> tuple[str, bool, GenotypeCall]:
    scenario, seed = task
    structures = [AlleleStructure.from_label(label) for label in scenario.alleles]
    alleles = [
        SimAllele(s, somatic_fraction=scenario.somatic if s.cag >= 36 else 0.0) for s in structures
    ]
    alleles[-1] = SimAllele(alleles[-1].structure, scenario.abundance, alleles[-1].somatic_fraction)
    sample = simulate(SimulationSpec(tuple(alleles), pairs=scenario.pairs, seed=seed))
    counts = count_reads(
        (a.sequence, b.sequence) for a, b in zip(sample.r1, sample.r2, strict=True)
    )
    call = call_genotype(counts)
    truth = true_genotype(structures)
    return scenario.name, call_matches([a.allele for a in call.alleles], truth), call


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, default=5)
    parser.add_argument("--workers", type=int)
    options = parser.parse_args()
    tasks = [(s, seed) for s in SCENARIOS for seed in range(options.seeds)]
    with ProcessPoolExecutor(options.workers) as pool:
        results = list(pool.map(run, tasks))

    print(f"{'scenario':<26}{'correct':>9}{'mean Q':>8}{'wrong, Q>=20':>14}  common flags")
    for scenario in SCENARIOS:
        rows = [(ok, call) for name, ok, call in results if name == scenario.name]
        correct = sum(ok for ok, _ in rows)
        confident_wrong = sum(not ok and call.quality >= 20 for ok, call in rows)
        mean_q = sum(call.quality for _, call in rows) / len(rows)
        flags = sorted({str(f) for _, call in rows for f in call.flags})
        print(
            f"{scenario.name:<26}{f'{correct}/{len(rows)}':>9}{mean_q:>8.1f}"
            f"{confident_wrong:>14}  {', '.join(flags)}"
        )
    for name, ok, call in results:
        if not ok:
            print(f"  wrong: {name}: called {call.label} (Q {call.quality:.1f})")


if __name__ == "__main__":
    main()
