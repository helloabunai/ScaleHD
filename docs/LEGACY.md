# Legacy genotyping (ScaleHD 1.x)

Pick "Legacy ScaleHD" as a job's genotyping method and each sample runs through
ScaleHD 1.x itself, completely unchanged. Version 1.3 as maintained by the Monckton lab
(MoncktonLab on GitHub, the copy in [`legacy/`](../legacy/), see `legacy/UPSTREAM`), on the
Python and libraries its `setup.py` pins. A sample gets what 1.x would give it when run on its
own, and `tools/legacy-parity.sh` checks that it does (see [DEVELOPMENT.md](DEVELOPMENT.md)).
That includes the lab's own test sample (`legacy/ContinuousTesting`), which gives exactly the
reports the lab expects.

1.x needs its own Python (3.12) and the programs it runs (bwa, samtools, cutadapt, FastQC,
seqtk, freebayes, picard), which is installed with the Docker image. If you run the new
stack outside of docker, it will only offer the model based genotype calling.
1.x also needs both read files of each sample, a job
won't run 1.x on a sample with only R1.

## Settings

The web interface provides the user to supply all configuration settings that would have been
manually typed into a XML config file with the old version.

You can set account default settings for legacy settings, if you wish.

## Cores per sample

The server has a number of cores for samples (`SCALEHD_WORKERS`, all of them by default).
A legacy job says how many each of its samples gets, which is the slider on the settings.
Where 1.x would give every sample every core on the machine, it gets this many (1.x's
`THREADS`, used by bwa, samtools, FastQC and its atypical allele scan).

So with e.g. 24 cores, 1 core per sample runs 24 samples at once, and 24 cores per sample runs
them one at a time (but faster). Fewer cores per job suits a batch of many samples,
more cores per job suits a batch of a few samples.

Sample core-queues are first come first served.

With the new model-based genotyping, one sample always uses one core (maybe up for debate)

## How a sample runs

Each sample from the new interface is run as a batch of one,
which is what 1.x does for each sample of a batch anyway, so a batch run and these runs 
give the same outputs. The sample's folder in the job's folder holds:

```
legacy-request.json   what the runner was asked to do for clarity. settings etc
legacy.log            1.x's log, and the output of the binaries it ran
legacy/
  config.xml          the configuration 1.x was given
  run/                1.x's own output: InstanceReport.csv, AlignedDistributions.csv,
                      ScaleHDLog.txt, the sample's SeqQC/Align/Predict folders, its graphs
  result.json         everything 1.x knew about the sample when it finished (used by the web UI)
```

1.x knows the sample by its name, or by its folder's name if its name has characters 1.x
can't use (anything but letters, digits, `_` and `-`). Demultiplexing renames it
`DMPX_<name>` in 1.x's own output, as 1.x always does.

1.x itself isn't changed. Around it, the runner (`apps/legacy-runner/scalehd_legacy_run.py`)
only sets `THREADS` to the cores per sample, skips the HTML report at the end (the web
interface shows the same and more), follows 1.x's log to report each stage as it starts,
and reads every value 1.x holds for the sample and its two alleles once it has finished.

## Results

- The genotype is 1.x's two alleles, primary (the one with most reads) first, as 1.x
  orders them.
- Confidence is 1.x's, a percentage for each allele: over 80% is good, over 60%
  believable, and less may need a look (1.x's documentation). The job's table shows the
  less confident allele's.
- Flags are those 1.x raised, by its InstanceReport column names. The sample's page
  lists all of them with what 1.x's documentation says each means, including "Heuristic
  Filtering Success", which says something went right and so isn't counted as a flag.
- The sample's page has each allele's values from 1.x's InstanceReport, the CAG and CCG
  read distributions 1.x charted, how the reads aligned, and the run (stages, time,
  cores, 1.x's version). Its files are there to download: everything 1.x knew about the
  sample (JSON), 1.x's PDF of graphs, and its log.
- A job's page has 1.x's InstanceReport.csv and AlignedDistributions.csv for the whole job,
  as 1.x writes them for a batch, with the samples in the job's order. The usual results
  export fills the columns 1.x has values for.
- When 1.x can't genotype a sample, the sample fails with 1.x's own reason and the stage
  it failed in, e.g. "Legacy ScaleHD failed in atypical scanning: <200 aligned reads in
  Allele #1. Data un-usable." 1.x's row is still in the job's InstanceReport, as 1.x
  writes it.

## Worth knowing

- When an allele's confidence is under 50, 1.x widens its confidence interval by a random
  1 or 2 CAG at the bottom and 1 at the top, without a fixed seed. Running the sample
  again can give a slightly different interval. The genotype doesn't change.
- 1.x runs some programs through a shell, so the workspace's path must not have spaces or
  other characters a shell treats specially. A sample in a workspace like that fails with
  a message saying so.
- Cancelling a job lets samples already running finish, as for any job.
