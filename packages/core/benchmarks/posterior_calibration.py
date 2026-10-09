"""Whether the caller's confidence can be believed (given limited data at time of writing).

A call with posterior 0.99 should be wrong about 1 time in 100. This simulates random
samples and compares each band of posteriors with how often those calls were right.
The samples are:

- genotypes as ``scalehd simulate-run`` draws them (60%), or neighbouring alleles one
  CAG apart (25%) or two apart (15%), the hard cases,
- at depths from 150 to 5,000 read pairs, so that some calls are unsure.

They are simulated once per pcr stutter preset (``SETTINGS``):

- ``calibrated``: the stutter the caller's priors come from,
- ``measured``: each sample's stutter drawn as far as samples in the ScaleHD 1.x training
  matrix vary (limited atm),
- ``prior``: drawn within the caller's prior spread, wider than that (what the priors
  themselves fit, so if the posteriors are calibrated anywhere it's here),
- ``twice``: double the prior spread (priors too narrow for the PCR).

Wrong calls at posterior 0.99 or more are listed with their flags.

Run with ``uv run python packages/core/benchmarks/posterior_calibration.py``.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np
from scalehd.counts import count_reads
from scalehd.genotype import CallerSettings, call_genotype
from scalehd.simulate import SimAllele, SimulationSpec, call_matches, simulate, true_genotype
from scalehd.simulate_run import placeholder_samples
from scalehd.structure import AlleleStructure
from stutter_variants import Variant

SETTINGS = {
    "calibrated": Variant("as calibrated"),
    "measured": Variant("drawn within the measured spread", measured=True),
    "prior": Variant("drawn within the prior spread", spread=1.0),
    "twice": Variant("drawn within double the prior spread", spread=2.0),
}

# Posterior bands: [low, high).
BANDS = ((0.0, 0.5), (0.5, 0.9), (0.9, 0.99), (0.99, 0.999), (0.999, 1.0 + 1e-9))


@dataclass(frozen=True)
class Sample:
    index: int
    alleles: tuple[str, ...]
    pairs: int


def draw_samples(n: int, seed: int) -> list[Sample]:
    rng = np.random.default_rng(seed)
    drawn = [s.alleles for s in placeholder_samples(n, seed) if s.name.startswith("random-")]
    samples = []
    for i in range(n):
        u = rng.random()
        if u < 0.6:
            alleles = drawn[i]
        else:
            cag = int(rng.integers(12, 46))
            gap = 1 if u < 0.85 else 2
            alleles = (f"{cag}_1_1_7_2", f"{cag + gap}_1_1_7_2")
        pairs = int(np.exp(rng.uniform(np.log(150), np.log(5000))))
        samples.append(Sample(i, alleles, pairs))
    return samples


@dataclass(frozen=True)
class Outcome:
    setting: str
    sample: Sample
    right: bool
    posterior: float
    called: str
    flags: tuple[str, ...]


def run(task: tuple[Sample, str, CallerSettings]) -> Outcome:
    sample, setting, settings = task
    structures = [AlleleStructure.from_label(label) for label in sample.alleles]
    stutter = SETTINGS[setting].curve(sample.index)
    spec = SimulationSpec(
        tuple(SimAllele(s) for s in structures),
        pairs=sample.pairs,
        stutter=stutter,
        seed=sample.index,
    )
    simulated = simulate(spec)
    counts = count_reads(
        (a.sequence, b.sequence) for a, b in zip(simulated.r1, simulated.r2, strict=True)
    )
    call = call_genotype(counts, settings)
    right = call_matches([a.allele for a in call.alleles], true_genotype(structures))
    flags = tuple(str(f) for f in call.flags)
    return Outcome(setting, sample, right, call.posterior, call.label, flags)


def report(setting: str, outcomes: list[Outcome]) -> None:
    rows = [(o.right, o.posterior) for o in outcomes]
    print(f"\nsimulated stutter {SETTINGS[setting].name}: {len(rows)} samples")
    print(f"  {'posterior':<14}{'samples':>9}{'mean posterior':>16}{'right':>9}")
    for low, high in BANDS:
        band = [(right, p) for right, p in rows if low <= p < high]
        if not band:
            print(f"  {f'{low:g}-{min(high, 1):g}':<14}{0:>9}")
            continue
        mean = sum(p for _, p in band) / len(band)
        share = sum(right for right, _ in band) / len(band)
        print(f"  {f'{low:g}-{min(high, 1):g}':<14}{len(band):>9}{mean:>16.4f}{share:>9.3f}")
    wrong = sum(not right for right, _ in rows)
    expected = sum(1 - p for _, p in rows)
    brier = sum((right - p) ** 2 for right, p in rows) / len(rows)
    print(f"  wrong {wrong}, expected from the posteriors {expected:.1f}, Brier score {brier:.4f}")
    for o in outcomes:
        if not o.right and o.posterior >= 0.99:
            print(
                f"  wrong at {o.posterior:.4f}: truth {'/'.join(o.sample.alleles):<24}"
                f"called {o.called:<24}{o.sample.pairs:>5} pairs  {','.join(o.flags)}"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=200, help="per stutter setting")
    parser.add_argument(
        "--settings", default=",".join(SETTINGS), help="stutter settings, comma-separated"
    )
    parser.add_argument("--local-radius", type=int, help="the caller's local_radius")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--workers", type=int)
    options = parser.parse_args()
    settings = CallerSettings()
    if options.local_radius is not None:
        settings = CallerSettings(local_radius=options.local_radius)
    chosen = options.settings.split(",")
    samples = draw_samples(options.samples, options.seed)
    tasks = [(sample, setting, settings) for setting in chosen for sample in samples]
    with ProcessPoolExecutor(options.workers) as pool:
        results = list(pool.map(run, tasks))
    for setting in chosen:
        report(setting, [o for o in results if o.setting == setting])


if __name__ == "__main__":
    main()
