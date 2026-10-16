"""ScaleHD 1.x genotyping. Not in this library: 1.x runs as itself, unchanged.

ScaleHD 1.x aligned reads with BWA-MEM to 4,000 synthetic references, re-aligned to
custom references when it found atypical alleles, and called genotypes from the
aligned read counts. This rework reads repeat structures straight from the reads and
has no alignment step. Rather than change code and bring 1.x's genotyper over, the server runs 1.x
itself on each sample in the Docker image (apps/legacy-runner, docs/LEGACY.md). See
MoncktonLab/ScaleHD on GitHub for the original (pinned in apps/legacy-runner/UPSTREAM).
"""

from __future__ import annotations

from typing import NoReturn


def call_genotype(*_: object, **__: object) -> NoReturn:
    raise NotImplementedError(
        "ScaleHD 1.x genotyping isn't in this library: the server runs 1.x itself in the"
        " Docker image (see docs/LEGACY.md)"
    )
