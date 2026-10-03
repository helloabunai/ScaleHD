"""ScaleHD 1.x genotyping. Not brought over yet.

ScaleHD 1.x aligned reads with BWA-MEM to 4,000 synthetic references, re-aligned to
custom references when it found atypical alleles, and called genotypes from the
aligned read counts. This rework reads repeat structures straight from the reads and
has no alignment step, so that has to come first. See ``legacy/`` at the repository
root for the original.
"""

from __future__ import annotations

from typing import NoReturn


def call_genotype(*_: object, **__: object) -> NoReturn:
    raise NotImplementedError("legacy (ScaleHD 1.x) genotyping is not available yet")
