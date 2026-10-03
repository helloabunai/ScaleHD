"""The demo job: simulated samples with known genotypes, through the real job machinery."""

from __future__ import annotations

from dataclasses import dataclass

from scalehd.simulate import true_genotype
from scalehd.structure import AlleleStructure

from .models import Job, Sample, User
from .schemas import GenotypeMethod, JobSettings


@dataclass(frozen=True)
class DemoSample:
    name: str
    alleles: tuple[str, ...]
    seed: int
    pairs: int = 5000


DEMO_SAMPLES: tuple[DemoSample, ...] = (
    DemoSample("normal-heterozygote", ("17_1_1_7_2", "21_1_1_7_2"), seed=1),
    DemoSample("expanded", ("17_1_1_7_2", "43_1_1_7_2"), seed=2),
    DemoSample("homozygous", ("21_1_1_7_2",), seed=3),
    DemoSample("neighbouring", ("17_1_1_7_2", "18_1_1_7_2"), seed=4),
    DemoSample("loss-of-interruption", ("42_0_1_7_2", "19_1_1_7_2"), seed=5),
    DemoSample("beyond-read-length", ("20_1_1_7_2", "95_1_1_7_2"), seed=6),
    DemoSample("ccg-7-and-10", ("17_1_1_7_2", "17_1_1_10_2"), seed=7),
    DemoSample("normal-ccg-10", ("19_1_1_10_2", "44_1_1_7_2"), seed=8),
    DemoSample("rare-ccg", ("16_1_1_6_2", "24_1_1_12_2"), seed=9),
    DemoSample("caacag-duplication", ("19_2_1_10_2", "40_1_1_7_2"), seed=10),
    DemoSample("ccgcca-deletion", ("19_1_0_7_2", "40_1_1_7_2"), seed=11),
    DemoSample("ccgcca-duplication", ("17_1_1_7_2", "42_1_2_7_2"), seed=12),
    DemoSample("ccgcca-deletion-and-insertion", ("19_1_0_7_2", "42_1_2_7_2"), seed=13),
)


def demo_job(owner: User, defaults: JobSettings) -> Job:
    """The demo job, not yet stored.

    Always the model method, because legacy can't run yet. The flag thresholds come
    from the user's defaults, so changing those changes the demo's flags.
    """
    settings = defaults.model_copy(update={"method": GenotypeMethod.MODEL, "call": True})
    job = Job(
        owner=owner,
        name=f"Demo: {len(DEMO_SAMPLES)} simulated samples",
        demo=True,
        settings=settings.model_dump(mode="json"),
    )
    job.samples = [
        Sample(
            name=sample.name,
            simulation={
                "alleles": list(sample.alleles),
                "pairs": sample.pairs,
                "seed": sample.seed,
            },
            truth=_truth(sample.alleles),
        )
        for sample in DEMO_SAMPLES
    ]
    return job


def _truth(alleles: tuple[str, ...]) -> str:
    structures = [AlleleStructure.from_label(label) for label in alleles]
    return "/".join(structure.label for structure in true_genotype(structures))
