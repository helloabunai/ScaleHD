# ScaleHD Rework

Genotyping of the Huntington disease *HTT* CAG/CCG repeat from paired-end amplicon
sequencing. This branch is a ground-up rewrite. The previous ScaleHD implementation
is kept for reference in [`legacy/`](legacy/). Looking back at old work is embarrassing
but everyone starts their professional journey somewhere. Hopefully this implementation
is more professional.

I am no longer working with the Monckton research group at University of Glasgow 
anymore so this is mostly just a hobby project that may or may not go anywhere.

## Quick start

Requires [uv](https://docs.astral.sh/uv/) and Python 3.13+.

```sh
uv sync
uv run scalehd simulate -a 17_1_1_7_2 -a 43_1_1_7_2 -n 20000 -o scratch --name s1
uv run scalehd count scratch/s1_R1.fastq.gz scratch/s1_R2.fastq.gz -o scratch/s1.counts.json
uv run scalehd call scratch/s1.counts.json -o scratch/s1.call.json
```

`scalehd genotype R1 R2` counts and calls in one step.

`simulate` writes paired FASTQ plus a `.truth.json`. Its model covers PCR stutter,
somatic expansion, length bias, sequencing errors that rise along the read, spacers
and adapter read-through, so the pipeline can be tested before real data is
available.

### Development

within `tools/` subdir is a mini script to refresh your docker dev stack. cos any dev is probably lazy

```sh
uv run pytest                        # tests
uv run ruff check packages apps      # lint
uv run ruff format packages apps     # format
uv run mypy                          # types
(cd apps/web && npm run build)       # frontend type-check and build
uv run python packages/core/benchmarks/parse_accuracy.py
uv run python packages/core/benchmarks/genotype_simulated.py
uv run python packages/core/benchmarks/legacy_matrix.py
```

`tools/check.sh` runs every check above in order (lint, format, types, frontend build, tests) and stops at the first failure. A git pre-push hook runs it before each push and
blocks the push if anything fails. Turn the hook in your local repo clone:

```sh
git config core.hooksPath tools/git-hooks
```

`git push --no-verify` skips it for one push.  Maybe don't do that.

#### From scratch: running the server with Docker

Everything (server, web interface, genotyping) runs in one container, so the machine
only needs Docker, Docker Compose and git.

1. Install Docker and Compose.

Example from scratch with ubuntu, except installing docker direct
from docker's repository as the OS repo is typically out of date.

```sh
sudo apt-get update
sudo apt-get install -y ca-certificates curl git
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc
echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] \
https://download.docker.com/linux/ubuntu $(. /etc/os-release && echo "$VERSION_CODENAME") stable" \
    | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
```

Then let your user run Docker without `sudo`, and log out and back in for it to take
effect (or run `newgrp docker` in the current shell):

```sh
sudo usermod -aG docker "$USER"
```

Membership of the `docker` group is root-equivalent on that machine, so only add
users who should have it.

Check both work:

```sh
docker run --rm hello-world
docker compose version
```

2. Get this code (currently on rework branch)

```sh
git clone https://github.com/helloabunai/scalehd.git
cd scalehd
git switch rework
```

3. Pick the folders and write `.env`.

The server needs two folders on the host. The container sees each at the same path as
the host does, so paths shown in the web interface are real host paths.

- A data folder with the sequencing runs (FASTQ files), mounted read-only. A network
  share works if it is mounted on the host before the container starts.
- A workspace folder for results, mounted read-write. ScaleHD will create a folder per user,
  and a subfolder per job. The root workspace/results folder must exist before starting. On Linux it must be writable by the container's user (the following paths are made up examples):

  ```sh
  sudo mkdir -p /srv/scalehd/workspace
  sudo chown 1000 /srv/scalehd/workspace
  ```

Specify such folder locations in the env file (copy the example template).

```sh
cp .env.example .env
```

Then edit `.env`. Only the first two settings are required at the moment:

| setting | default | meaning |
|---|---|---|
| `SCALEHD_DATA_ROOT` | required | sequencing data users can browse, read-only (absolute path) |
| `SCALEHD_WORKSPACE` | required | where results go (absolute path, must exist) |
| `SCALEHD_WORKERS` | every core | samples processed at once, one process each |
| `SCALEHD_ALLOW_REGISTRATION` | `true` | whether anyone who can reach the server may create an account (the first account, the admin, can always be created) |
| `SCALEHD_SESSION_DAYS` | `30` | how long a login lasts |

The database lives in a Docker volume (`scalehd_scalehd-data`), separate from the
workspace. `SCALEHD_DATABASE_DIR` and `SCALEHD_WEB_DIR` are set by the image; leave
them alone.

Most settings are work in progress:

- The data folder is mounted but can't be browsed from the web interface yet, so
  only the demo job (simulated samples) runs.
- Of the genotyping methods (Settings page, and per job), Legacy (ScaleHD 1.x) is not
  available yet and New (model-based) is a beta. The demo always uses New.
- The flag thresholds and the other job settings are not fully implemented yet so ignore those.

4. Start it

```sh
docker compose up -d --build
docker compose ps        # STATUS should become "healthy"
docker compose logs -f   # server log (Ctrl+C stops following, not the server)
```

The container restarts on its own after a reboot (`restart: unless-stopped`).
`docker compose down` stops it; the database volume and the folders are kept.

5. Open the web interface.

On the same machine: <http://localhost:8000>. The first account you register is
the admin. From there, "Run demo" on the home page runs the simulated samples
end to end. If you were running the docker on a lab computer then
you'd access it at <http://whatever_dedicated_ip_address:8000> from the same internal network. 

Network and port customisation will be implemented at some point in the future.

### Updating

```sh
git pull
docker compose up -d --build
```

There are no database migrations yet. If the server refuses to start because the
database is from an older version (`docker compose logs` says so), `tools/refresh.sh`
deletes the database volume and rebuilds. All accounts and job history are lost;
files in the workspace are kept. Development work in progress baby !!!

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
confidence, was rough in the previous implementation. Now each candidate genotype is a
mixture of two alleles, each smeared by a PCR stutter kernel whose shape depends on
repeat length (calibrated on the ScaleHD 1.x training matrix). Candidates are compared
by how well they explain every molecule, which gives a posterior probability for the
call plus flags for the cases that deserve a look. See `packages/core/src/scalehd/genotype.py`.

## Layout

```
packages/core/        scalehd: pure-Python library and CLI (no web or DB dependencies)
  src/scalehd/        structure, amplicon, parse, pairs, counts, calibration, genotype,
                      simulate, seqio, cli
  tests/              unit, property-based and end-to-end tests
  benchmarks/         accuracy on simulated data and the legacy labelled matrix
apps/server/          scalehd-server: FastAPI, job runner, SQLite database (skeleton)
apps/web/             web interface: React, TypeScript, Vite (skeleton)
Dockerfile            one image with the server, the core and the built frontend
compose.yaml          runs that image with a data volume and a read-only FASTQ folder
legacy/               ScaleHD 1.x, for reference only
```

## Using real FASTQ files

### What the input should look like

Basically the same rules as the previous implementation of ScaleHD because i know nothing else.

- *Paired-end amplicon reads, one R1 and one R2 file per sample, gzipped or not.
  Single-end works too (pass R1 only), but you obviously lose R2 cross-checking and the
  CCG/CCT side of long alleles.
- R1 should be for the CAG strand, starting at the 5' flank and runs into the CAG tract,
  as with the GeM-HD MiSeq primers. R2 runs the other way and is reverse-complement.
- Read pairs in step. Record *n* of R1 must be record *n* of R2. Files
  straight off the sequencer are fine (cos that's what we had). If you filter or trim 
  R1 and R2 separately and drop reads from one file only, pairs fall out of step and
  `scalehd` will stop with an error.
- No trimming or demultiplexing needed (different to 1.0). Mostly borne out of the fact
  I am simulating data for now so this may change. Reads are located by 20-base anchors either side of 
  the repeat (`TCGAGTCCCTCAAGTCCTTC` 5', `CAGCTTCCTCAGCCGCCGCC` 3'),
  so primers, heterogeneity spacers, adapters and reads from other amplicons in the
  same library are skipped. If your primers sit inside those anchors, the default
  amplicon won't fit (see below).
- All work in progress and let's see if I even finish it. If I get access to real data again
  then this may all change.

### One sample

```sh
uv run scalehd genotype sample_R1.fastq.gz sample_R2.fastq.gz \
    --counts sample.counts.json -o sample.call.json
```

This prints the call and writes two files:

- `sample.counts.json`: every molecule's repeat structure, tallied. Calling can be
  rerun from this alone with `uv run scalehd call sample.counts.json`.
- `sample.call.json`: the genotype, each allele's structure, share of molecules,
  fitted stutter, the ScaleHD 1.x slippage and mosaicism ratios, flags, and the
  runner-up genotypes with their probabilities.

### A batch of samples

Assuming files named `<sample>_R1.fastq.gz` / `<sample>_R2.fastq.gz`:

```sh
mkdir -p results
for r1 in run/*_R1.fastq.gz; do
    sample=$(basename "$r1" _R1.fastq.gz)
    uv run scalehd genotype "$r1" "run/${sample}_R2.fastq.gz" \
        --counts "results/$sample.counts.json" -o "results/$sample.call.json" \
        > "results/$sample.txt" &
done
wait
```

TODO: improve input to have folder support to make this easier to use.

Each sample is independent, so running them in parallel (the trailing `&`) is safe.
On a machine with many cores, cap the number running at once (for example with
`xargs -P`) rather than starting hundreds together.

### Reading the result

- Posterior / quality: the probability that the genotype is right given the
  model, and the same thing on a Phred scale (quality 20 = 1 in 100 wrong, capped
  at 99). It covers stutter and sampling noise. It does not cover the stutter model
  itself being wrong for PCR conditions (which again.. simulated data lol).
- Allele labels are same as ScaleHD 1.x `CAG_CAACAG_CCGCCA_CCG_CCT`.
  `42_0_1_7_2` is a loss of the CAA interruption, `19_2_1_10_2` a CAACAG duplication.
  `83+_1_1_7_2` means no read spanned the CAG tract, so only a lower bound is known.
  A rough estimate of the true length is given only when the reads also limit it from
  above, and it leans on the stutter model beyond the lengths it was measured at. On
  simulated data the reads never do, so expect the lower bound alone.
- Flags mark what deserves a manual inspection (unchanged really from previous):

  | flag | meaning |
  |---|---|
  | `low_depth` | fewer than 500 usable molecules |
  | `low_confidence` | posterior below 0.99 |
  | `homozygous` | both alleles identical |
  | `neighbouring` | alleles one CAG apart. the hardest case to separate from stutter |
  | `atypical` | an allele without the common `1_1_x_2` intervening and CCT structure |
  | `beyond_read_length` | an allele longer than the reads. CAG is a lower bound not definitive |
  | `allele_imbalance` | one allele has under 20% of molecules |
  | `high_background` | over 5% of molecules fit neither allele |
  | `unexplained_peak` | a peak the genotype doesn't explain. third allele, contamination, mosaicism? |
  | `high_discordance` | over 15% of molecules dropped because their read base pairings disagreed |

### Checking the input went well

The first lines `scalehd genotype` prints (and `read_outcomes` in the counts JSON)
show how reads fared:

- mostly `no_anchor` in R1 = the files are probably swapped (R1 is the CCG strand)
  or the amplicon doesn't contain the default anchors.
- many `nonconforming` = reads reach both anchors but don't fit the tract structure.
  Expect this with very poor quality or human error
- high `dropped` = read pairs often disagree. A few percent is normal sequencing error,
  much more suggests quality problems or out-of-step bases.

### Other amplicons and PCR conditions

- Other flanks. The default flanks are those of the ScaleHD 1.x reference library.
  From Python, `RepeatParser(amplicon=AmpliconSpec(...))` and `count_fastq(...,
  parser=...)` take other flanks. The CLI doesn't expose this yet i.e. TODO work
- Other stutter Stutter priors were measured on MiSeq data from one PCR protocol
  (`scalehd.calibration.HTT_MISEQ`). Very different chemistry or cycle numbers may need
  a recalibrated curve, passed via `CallerSettings(stutter=...)`. also still WIP.

## Web interface and API (skeleton)

Basic interface and user account system. Everything else answers "not implemented yet".
See [`apps/server/README.md`](apps/server/README.md) for what is real and what's a stub.

Don't bother running anything. Massively WIP.

## Licence

MIT, see [LICENSE](LICENSE).
