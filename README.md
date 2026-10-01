# ScaleHD Rework

Genotyping of the Huntington disease *HTT* CAG/CCG repeat from paired-end amplicon
sequencing. This branch is a ground-up rewrite. The previous ScaleHD implementation
is kept for reference in [`legacy/`](legacy/). Looking back at old work is embarrassing
but everyone starts their professional journey somewhere. Hopefully this implementation
is more professional.

I am no longer working with the Monckton research group at University of Glasgow 
anymore so this is mostly just a hobby project that may or may not go anywhere.

## What's different

The original implementation of ScaleHD was a command line python package, which took settings
from users via an XML file, and executed pipeline steps based on what was requested.
It would align sequencing reads with BWA-MEM to 4,000 synthetic references,
scanned the presented structure for atypical alleles, and re-aligned to custom references when
any were found. Genotyping was done with fuzzy logic that in hindsight could be improved 
a lot.

I aim for this re-work to be a docker container (for ease of shipping), which provides a 
web based front-end, where users can specify settings, view past jobs/results, run specific
ScaleHD features based on requirements, and export results to PDF files. While there will be an
API for calling backend functionality from the web interface, the backend/python package
can also be used as before, i.e. a command line interface, if users prefer.

I'm also messing around with the logic for calling the genotypes which means results may be 
compeltely incomparable but I feel that alignment is perhaps not required so much for this
task.

ScaleHD rework, or ScaleHD2, or whatever this is called, reads the repeat structure straight 
from each read instead:

1. Short sequences at the flank ends next to the repeat locate it in the
   read. Primers, spacers, adapters and off-target reads should thus need no trimming first.
2. The bases between the anchors are split into the established Huntingon structure
   `(CAG)n (CAACAG)a (CCGCCA)b (CCG)m (CCT)k`, tolerating substitutions and single-base
   insertions or deletions. 
3. R1 and R2 are the same molecule. When both see a repeat tract and
   disagree, the molecule is dropped, because sequencing errors rarely coincide while
   PCR stutter is shared. For long alleles, R1 supplies the CAG tract and R2 the CCG
   and CCT tracts.
4. A read that ends inside the repeat only gives a lower bound (`83+`).
   It is never forced onto a reference.

One major issue with me attempting to do this re-write is that because I'm no longer at the
university, I have no access to test data. But I'm bored so I'm doing it anyway. I've written a
simulator to simulate sequencing reads based on my unreliable memory of Huntington disease.

Nobody should probably quote any of the science that I remember/have read from our papers.

Genotype calling, the stage that turns per-molecule counts into two alleles with a
confidence, os not implemented yet. But the previous implemention was, eh.. rough..

## Layout

```
packages/core/        scalehd: pure-Python library and CLI (no web or DB dependencies)
  src/scalehd/        structure, amplicon, parse, pairs, counts, simulate, seqio, cli
  tests/              unit, property-based and end-to-end tests
  benchmarks/         accuracy and speed on simulated data
legacy/               ScaleHD 1.x, for reference only
```

Planned: `apps/server` (FastAPI, jobs, database) and `apps/web` (frontend).

## Quick start

Requires [uv](https://docs.astral.sh/uv/) and Python 3.13+.

```sh
uv sync
uv run scalehd simulate -a 17_1_1_7_2 -a 43_1_1_7_2 -n 20000 -o scratch --name s1
uv run scalehd count scratch/s1_R1.fastq.gz scratch/s1_R2.fastq.gz -o scratch/s1.counts.json
```

`simulate` writes paired FASTQ plus a `.truth.json`. Its model covers PCR stutter,
somatic expansion, length bias, sequencing errors that rise along the read, spacers
and adapter read-through, so the pipeline can be tested before real data is
available.

## Development

```sh
uv run pytest                 # tests
uv run ruff check packages    # lint
uv run ruff format packages   # format
uv run mypy                   # types
uv run python packages/core/benchmarks/parse_accuracy.py
```

## Licence

MIT, see [LICENSE](LICENSE).
