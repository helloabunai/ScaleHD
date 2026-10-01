"""Simulate paired-end HTT amplicon reads from known genotypes.

The model is simple but covers the artefacts a caller has to cope with (hopefully):

- PCR stutter, mostly contractions and growing with repeat length (calibrated on the
  ScaleHD 1.x training matrix, see `scalehd.calibration`)
- somatic expansion
- amplification bias toward shorter alleles
- sequencing errors that rise along the read
- heterogeneity spacers
- adapter read-through

I may have forgotten many science words !! yipee !!

Such simulations are completely made up by me, the idiot, so probably not "good".
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from .amplicon import HTT_AMPLICON, AmpliconSpec
from .calibration import HTT_MISEQ, StutterCurve, cag_kernel
from .seqio import FastqRecord, open_text, reverse_complement, write_fastq
from .structure import AlleleStructure

SCHEMA = "scalehd.simulation/1"

NEXTERA_READTHROUGH_R1 = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
NEXTERA_READTHROUGH_R2 = "CTGTCTCTTATACACATCTGACGCTGCCGACGA"

_BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
_BASE_INDEX = np.zeros(256, dtype=np.uint8)
_BASE_INDEX[_BASES] = np.arange(4, dtype=np.uint8)


@dataclass(frozen=True, slots=True)
class SimAllele:
    structure: AlleleStructure
    abundance: float = 1.0
    # Fraction of templates carrying a somatic expansion, and the mean extra CAG units
    # (geometrically distributed, at least one) among those that do.
    somatic_fraction: float = 0.0
    somatic_mean: float = 3.0

    def __post_init__(self) -> None:
        if self.abundance <= 0:
            raise ValueError("abundance must be positive")
        if not 0 <= self.somatic_fraction <= 1:
            raise ValueError("somatic_fraction must be between 0 and 1")
        if self.somatic_mean < 1:
            raise ValueError("somatic_mean must be at least 1")


@dataclass(frozen=True, slots=True)
class SequencingModel:
    read_length: int = 300
    # Substitution rate rises from error_start to error_end along the read.
    error_start: float = 0.001
    error_end: float = 0.02
    error_shape: float = 3.0
    max_spacer: int = 3
    readthrough_r1: str = NEXTERA_READTHROUGH_R1
    readthrough_r2: str = NEXTERA_READTHROUGH_R2

    def error_profile(self) -> np.ndarray:
        position = np.linspace(0.0, 1.0, self.read_length)
        return self.error_start + (self.error_end - self.error_start) * position**self.error_shape

    def quality_string(self) -> str:
        phred = np.clip(np.round(-10 * np.log10(self.error_profile())), 2, 41).astype(int)
        return "".join(chr(33 + q) for q in phred)


@dataclass(frozen=True, slots=True)
class SimulationSpec:
    alleles: tuple[SimAllele, ...]
    pairs: int = 20_000
    # Relative amplification penalty per CAG unit above the shortest allele. Real samples
    # show little i.e. the shorter allele's median share of molecules is 0.48.
    length_bias: float = 0.0
    off_target: float = 0.0
    stutter: StutterCurve = HTT_MISEQ
    # Probability of each of a one-unit CCG contraction and expansion.
    ccg_slippage: float = 0.01
    sequencing: SequencingModel = field(default_factory=SequencingModel)
    amplicon: AmpliconSpec = HTT_AMPLICON
    seed: int = 0

    def __post_init__(self) -> None:
        if not self.alleles:
            raise ValueError("at least one allele is required")


@dataclass
class SimulatedSample:
    spec: SimulationSpec
    r1: list[FastqRecord]
    r2: list[FastqRecord]
    molecules: Counter[AlleleStructure]
    off_target_pairs: int

    def truth(self, name: str) -> dict[str, Any]:
        spec = self.spec
        return {
            "schema": SCHEMA,
            "sample": name,
            "seed": spec.seed,
            "pairs": spec.pairs,
            "alleles": [
                {
                    "structure": a.structure.label,
                    "abundance": a.abundance,
                    "somatic_fraction": a.somatic_fraction,
                    "somatic_mean": a.somatic_mean,
                }
                for a in spec.alleles
            ],
            "length_bias": spec.length_bias,
            "off_target": spec.off_target,
            "stutter": asdict(spec.stutter),
            "ccg_slippage": spec.ccg_slippage,
            "sequencing": asdict(spec.sequencing),
            "amplicon": asdict(spec.amplicon),
            "molecules": [
                {"structure": s.label, "count": n} for s, n in self.molecules.most_common()
            ],
            "off_target_pairs": self.off_target_pairs,
        }

    def write(self, directory: str | Path, name: str) -> tuple[Path, Path, Path]:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        r1_path = directory / f"{name}_R1.fastq.gz"
        r2_path = directory / f"{name}_R2.fastq.gz"
        truth_path = directory / f"{name}.truth.json"
        for path, records in ((r1_path, self.r1), (r2_path, self.r2)):
            with open_text(path, "wt") as handle:
                write_fastq(handle, records)
        truth_path.write_text(json.dumps(self.truth(name), indent=2) + "\n")
        return r1_path, r2_path, truth_path


def simulate(spec: SimulationSpec) -> SimulatedSample:
    rng = np.random.default_rng(spec.seed)
    seq = spec.sequencing
    profile = seq.error_profile()
    quality = seq.quality_string()

    shortest = min(a.structure.cag for a in spec.alleles)
    weights = np.array(
        [a.abundance * (1 - spec.length_bias) ** (a.structure.cag - shortest) for a in spec.alleles]
    )
    off_target_pairs = int(rng.binomial(spec.pairs, spec.off_target))
    template = rng.choice(
        len(spec.alleles), size=spec.pairs - off_target_pairs, p=weights / weights.sum()
    )

    molecules: Counter[AlleleStructure] = Counter()
    r1: list[FastqRecord] = []
    r2: list[FastqRecord] = []

    def read(insert: str, readthrough: str) -> str:
        spacer = _random_bases(rng, int(rng.integers(0, seq.max_spacer + 1)))
        body = spacer + insert + readthrough
        if len(body) < seq.read_length:
            body += _random_bases(rng, seq.read_length - len(body))
        return _add_errors(rng, body[: seq.read_length], profile)

    for i, index in enumerate(template):
        allele = spec.alleles[index]
        molecule = _amplify(rng, allele, spec.stutter, spec.ccg_slippage)
        molecules[molecule] += 1
        insert = spec.amplicon.sequence(molecule.repeat_sequence())
        r1.append(
            FastqRecord(f"sim{i} 1:N:0:{molecule.label}", read(insert, seq.readthrough_r1), quality)
        )
        r2.append(
            FastqRecord(
                f"sim{i} 2:N:0:{molecule.label}",
                read(reverse_complement(insert), seq.readthrough_r2),
                quality,
            )
        )

    for i in range(len(template), spec.pairs):
        for records, mate in ((r1, 1), (r2, 2)):
            junk = _random_bases(rng, seq.read_length)
            records.append(FastqRecord(f"sim{i} {mate}:N:0:off_target", junk, quality))

    return SimulatedSample(spec, r1, r2, molecules, off_target_pairs)


def _amplify(
    rng: np.random.Generator, allele: SimAllele, stutter: StutterCurve, ccg_slippage: float
) -> AlleleStructure:
    """One PCR product. optional extra somatic expansion, then stutter on template."""
    structure = allele.structure
    cag = structure.cag
    if allele.somatic_fraction and rng.random() < allele.somatic_fraction:
        cag += int(rng.geometric(1 / allele.somatic_mean))
    shifts, probabilities = cag_kernel(stutter, cag)
    cag = max(1, cag + int(rng.choice(shifts, p=probabilities)))
    ccg = structure.ccg
    if (u := rng.random()) < ccg_slippage:
        ccg -= 1
    elif u < 2 * ccg_slippage:
        ccg += 1
    return structure.with_counts(cag=cag, ccg=max(1, ccg))


def _random_bases(rng: np.random.Generator, n: int) -> str:
    return _BASES[rng.integers(0, 4, size=n)].tobytes().decode("ascii")


def _add_errors(rng: np.random.Generator, sequence: str, profile: np.ndarray) -> str:
    bases = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8).copy()
    hits = np.flatnonzero(rng.random(bases.size) < profile[: bases.size])
    if hits.size:
        shifted = (_BASE_INDEX[bases[hits]] + rng.integers(1, 4, size=hits.size)) % 4
        bases[hits] = _BASES[shifted]
    return bases.tobytes().decode("ascii")
