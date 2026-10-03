"""Genotype calling: per-molecule repeat structures in, a genotype call out.

Each method is a subpackage:

- ``model``: the model-based caller written for this rework (see its ``__init__``).
- ``legacy``: ScaleHD 1.x genotyping. Not brought over yet.

The model method's API is re-exported here, so ``from scalehd.genotype import
call_genotype`` is the model-based caller.
"""

from . import legacy, model
from .model import (
    SCHEMA,
    AlleleCall,
    CallerSettings,
    Candidate,
    Flag,
    GenotypeCall,
    NoMoleculesError,
    call_genotype,
    candidate_alleles,
)

__all__ = [
    "SCHEMA",
    "AlleleCall",
    "CallerSettings",
    "Candidate",
    "Flag",
    "GenotypeCall",
    "NoMoleculesError",
    "call_genotype",
    "candidate_alleles",
    "legacy",
    "model",
]
