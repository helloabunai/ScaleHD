"""Genotype calls on simulated samples whose stutter isn't what the caller was modelled on.

``genotype_simulated.py`` simulates with the same stutter curve the caller's priors come
from, so it can't show what happens when real PCR stutters differently. Here every
scenario of that benchmark is simulated with each stutter variant of
``stutter_variants.py`` (less or more stutter, other tails, the curve moved along the CAG
axis, random per sample), and called with the caller's usual priors.

For each variant, how many calls are right, how many are wrong at quality 20 or more
(the ones nobody would check maybe), and how many wrong calls the posteriors expected (the
sum of 1 - posterior). Each wrong call is listed with its flags. Run with
``uv run python packages/core/benchmarks/stutter_robustness.py``.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

from genotype_simulated import SCENARIOS, Scenario
from scalehd.counts import count_reads
from scalehd.genotype import call_genotype
from scalehd.simulate import SimAllele, SimulationSpec, call_matches, simulate, true_genotype
from scalehd.structure import AlleleStructure
from stutter_variants import VARIANTS, Variant


@dataclass(frozen=True)
class Outcome:
    variant: str
    scenario: str
    seed: int
    right: bool
    called: str
    posterior: float
    quality: float
    flags: tuple[str, ...]


def run(task: tuple[Variant, Scenario, int]) -> Outcome:
    variant, scenario, seed = task
    # A random variant draws its own curve for every sample, not one per seed.
    draw = seed * len(SCENARIOS) + SCENARIOS.index(scenario)
    structures = [AlleleStructure.from_label(label) for label in scenario.alleles]
    alleles = [
        SimAllele(s, somatic_fraction=scenario.somatic if s.cag >= 36 else 0.0) for s in structures
    ]
    alleles[-1] = SimAllele(alleles[-1].structure, scenario.abundance, alleles[-1].somatic_fraction)
    spec = SimulationSpec(
        tuple(alleles), pairs=scenario.pairs, stutter=variant.curve(draw), seed=seed
    )
    sample = simulate(spec)
    counts = count_reads(
        (a.sequence, b.sequence) for a, b in zip(sample.r1, sample.r2, strict=True)
    )
    call = call_genotype(counts)
    right = call_matches([a.allele for a in call.alleles], true_genotype(structures))
    flags = tuple(str(f) for f in call.flags)
    return Outcome(
        variant.name, scenario.name, seed, right, call.label, call.posterior, call.quality, flags
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, default=3)
    parser.add_argument("--workers", type=int)
    options = parser.parse_args()
    tasks = [(v, s, seed) for v in VARIANTS for s in SCENARIOS for seed in range(options.seeds)]
    with ProcessPoolExecutor(options.workers) as pool:
        outcomes = list(pool.map(run, tasks))

    print(f"{'simulated stutter':<28}{'right':>10}{'wrong, Q>=20':>14}{'expected wrong':>16}")
    for variant in VARIANTS:
        rows = [o for o in outcomes if o.variant == variant.name]
        right = sum(o.right for o in rows)
        confident = sum(not o.right and o.quality >= 20 for o in rows)
        expected = sum(1 - o.posterior for o in rows)
        print(f"{variant.name:<28}{f'{right}/{len(rows)}':>10}{confident:>14}{expected:>16.1f}")

    print("\nwrong calls:")
    for o in outcomes:
        if not o.right:
            truth = next(s for s in SCENARIOS if s.name == o.scenario).alleles
            print(
                f"  {o.variant:<28}{o.scenario:<26}seed {o.seed}  "
                f"truth {'/'.join(truth):<24}called {o.called:<24}Q {o.quality:5.1f}  "
                f"{','.join(o.flags)}"
            )


if __name__ == "__main__":
    main()
