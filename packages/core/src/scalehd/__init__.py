"""ScaleHD: HTT CAG/CCG repeat genotyping from amplicon sequencing."""

import os

# Single-threaded by default, unless the user has chosen otherwise.
for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_variable, "1")

from importlib.metadata import version

from .amplicon import HTT_AMPLICON, AmpliconSpec
from .counts import SampleCounts, count_fastq, count_reads
from .parse import ParserSettings, ReadOutcome, ReadParse, RepeatParser
from .structure import AlleleStructure, FieldStatus, Observation

__version__ = version("scalehd")

__all__ = [
    "HTT_AMPLICON",
    "AlleleStructure",
    "AmpliconSpec",
    "FieldStatus",
    "Observation",
    "ParserSettings",
    "ReadOutcome",
    "ReadParse",
    "RepeatParser",
    "SampleCounts",
    "__version__",
    "count_fastq",
    "count_reads",
]
