"""SNP calls on simulated samples (single-base substitutions only for now).

Each scenario is simulated with GeM-HD primers and several seeds, with PCR copying errors
as well as sequencing errors, then counted and called from its reads. A SNP is found when
one is called at the true site with the true base. It is on the right allele when the
call puts it on the allele that carries it (alleles in order, shorter first) i.e. "allele
unknown" when the call couldn't tell the alleles apart.
Run with ``uv run python packages/core/benchmarks/snp_simulated.py``.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

from scalehd.amplicon import HTT_GEM_HD
from scalehd.counts import count_reads
from scalehd.genotype import GenotypeCall, call_genotype
from scalehd.parse import RepeatParser
from scalehd.simulate import SimAllele, SimulationSpec, simulate
from scalehd.structure import AlleleStructure
from scalehd.substitutions import FIVE_PRIME, Site

type Snps = tuple[tuple[Site, str], ...]


@dataclass(frozen=True)
class Scenario:
    name: str
    # Each allele's structure and the SNPs it carries.
    alleles: tuple[tuple[str, Snps], ...]
    pairs: int = 5000
    pcr_error: float = 1e-4


def _cag(unit: int, offset: int, base: str) -> Snps:
    return ((Site("cag", unit, offset), base),)


def _ccg(unit: int, offset: int, base: str) -> Snps:
    return ((Site("ccg", unit, offset), base),)


FLANK = ((Site(FIVE_PRIME, 3), "A"),)
SCENARIOS = (
    Scenario("5' flank, short allele", (("17_1_1_7_2", FLANK), ("43_1_1_7_2", ()))),
    Scenario("CAG, expanded allele", (("17_1_1_7_2", ()), ("43_1_1_7_2", _cag(12, 2, "T")))),
    Scenario("CAG near its end", (("17_1_1_7_2", ()), ("43_1_1_7_2", _cag(41, 2, "T")))),
    Scenario("CCG, short allele", (("17_1_1_7_2", _ccg(3, 0, "T")), ("21_1_1_7_2", ()))),
    Scenario("one on each allele", (("17_1_1_7_2", FLANK), ("43_1_1_7_2", _ccg(4, 0, "T")))),
    Scenario("same length, one copy", (("20_1_1_7_2", ()), ("20_1_1_7_2", _ccg(3, 0, "T")))),
    Scenario("same length, both copies", (("20_1_1_7_2", _cag(5, 2, "A")),) * 2),
    Scenario("neighbouring", (("17_1_1_7_2", ()), ("18_1_1_7_2", _cag(5, 2, "A")))),
    Scenario("long, 66", (("20_1_1_7_2", ()), ("66_1_1_7_2", _cag(30, 2, "T")))),
    Scenario("very long, 75, CCG", (("20_1_1_7_2", ()), ("75_1_1_7_2", _ccg(4, 2, "A")))),
    Scenario("beyond read length, CCG", (("20_1_1_7_2", ()), ("95_1_1_7_2", _ccg(4, 2, "A")))),
    Scenario("low depth", (("17_1_1_7_2", ()), ("43_1_1_7_2", _cag(12, 2, "T"))), pairs=300),
    Scenario("none", (("17_1_1_7_2", ()), ("43_1_1_7_2", ()))),
    Scenario("none, long", (("20_1_1_7_2", ()), ("70_1_1_7_2", ()))),
    Scenario("none, more PCR errors", (("17_1_1_7_2", ()), ("43_1_1_7_2", ())), pcr_error=5e-4),
)


@dataclass(frozen=True)
class Outcome:
    expected: int
    found: int
    right_allele: int
    unknown_allele: int
    false: int
    qualities: tuple[float, ...]
    call: GenotypeCall


def run(task: tuple[Scenario, int]) -> tuple[str, Outcome]:
    scenario, seed = task
    alleles = tuple(
        SimAllele(AlleleStructure.from_label(label), snps=snps) for label, snps in scenario.alleles
    )
    spec = SimulationSpec(
        alleles, pairs=scenario.pairs, seed=seed, amplicon=HTT_GEM_HD, pcr_error=scenario.pcr_error
    )
    sample = simulate(spec)
    counts = count_reads(zip(sample.r1, sample.r2, strict=True), RepeatParser(HTT_GEM_HD))
    call = call_genotype(counts)

    # Which alleles carry each true SNP, by their place in the scenario (shorter first, as
    # in the call).
    carriers: dict[tuple[Site, str], set[int]] = {}
    for index, (_, snps) in enumerate(scenario.alleles):
        for snp in snps:
            carriers.setdefault(snp, set()).add(index)
    found = right = unknown = false = 0
    qualities = []
    for snp in call.snps or ():
        key = (snp.site, snp.alt)
        if key not in carriers:
            false += 1
            continue
        found += 1
        qualities.append(snp.quality)
        if snp.allele is None:
            unknown += 1
        elif snp.allele in carriers[key]:
            right += 1
    outcome = Outcome(len(carriers), found, right, unknown, false, tuple(qualities), call)
    return scenario.name, outcome


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, default=3)
    parser.add_argument("--workers", type=int)
    options = parser.parse_args()
    tasks = [(s, seed) for s in SCENARIOS for seed in range(options.seeds)]
    with ProcessPoolExecutor(options.workers) as pool:
        results = list(pool.map(run, tasks))

    print(
        f"{'scenario':<28}{'found':>8}{'right allele':>14}{'allele unknown':>16}"
        f"{'false':>7}{'mean Q':>8}"
    )
    for scenario in SCENARIOS:
        rows = [outcome for name, outcome in results if name == scenario.name]
        expected = sum(o.expected for o in rows)
        found = sum(o.found for o in rows)
        qualities = [q for o in rows for q in o.qualities]
        mean_q = f"{sum(qualities) / len(qualities):.1f}" if qualities else "-"
        print(
            f"{scenario.name:<28}{f'{found}/{expected}':>8}{sum(o.right_allele for o in rows):>14}"
            f"{sum(o.unknown_allele for o in rows):>16}{sum(o.false for o in rows):>7}{mean_q:>8}"
        )
    for name, outcome in results:
        if outcome.found < outcome.expected or outcome.false:
            called = "; ".join(
                f"{snp.label} on allele {snp.allele} (Q {snp.quality:.0f})"
                for snp in outcome.call.snps or ()
            )
            print(f"  {name}: genotype {outcome.call.label}, SNPs: {called or 'none'}")


if __name__ == "__main__":
    main()
