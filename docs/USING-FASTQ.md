# Using real FASTQ files

## What the input should look like

Basically the same rules as the previous implementation of ScaleHD because i know nothing else.

- Paired-end amplicon reads, one R1 and one R2 file per sample, gzipped or not.
  Single-end works too (pass R1 only), but you obviously lose R2 cross-checking and the
  CCG/CCT side of long alleles, which will affect genotyping results.
- R1 should be primarily for the CAG strand, starting at the 5' flank and into the CAG tract,
  as with the GeM-HD MiSeq primers. R2 runs the other way and is reverse-complement.
- Read pairs in step. Record *n* of R1 must be record *n* of R2. Files
  straight off the sequencer are fine (cos that's what we had???). If you filter or trim
  R1 and R2 separately and drop reads from one file only, pairs fall out of step and
  `scalehd` will stop with an error.
- No trimming or demultiplexing needed (different to 1.0). Mostly borne out of the fact
  I am simulating data for now so this may change. Reads are located by 20-base anchors either side of
  the repeat (`TCGAGTCCCTCAAGTCCTTC` 5', `CAGCTTCCTCAGCCGCCGCC` 3'),
  so primers, heterogeneity spacers, adapters and reads from other amplicons in the
  same library are skipped. If your primers sit inside those anchors, the default
  amplicon won't fit (see below).
- All work in progress and subject to change. Feedback from real scientists will improve my assumptions
  and may alter how certain features work.

## Processing one sample

```sh
uv run scalehd genotype sample_R1.fastq.gz sample_R2.fastq.gz \
    --counts sample.counts.json -o sample.call.json
```

This prints the call and writes two files:

- `sample.counts.json`: every molecule's repeat structure, tallied. Calling can be
  rerun from this alone with `uv run scalehd call sample.counts.json`.
- `sample.call.json`: the automated genotype call, each allele's structure, share of molecules,
  fitted PCR stutter, the ScaleHD 1.x slippage and mosaicism ratios, analysis flags, and the
  runner-up genotypes with their probabilities for manual inspection.

## Processing a batch of samples

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

## Reading the result

- Posterior / quality: the probability that the genotype is right given the
  model, and the same thing on a Phred scale (quality 20 = 1 in 100 wrong, capped
  at 99). It covers PCR stutter and sampling noise. It does not cover the stutter model
  itself being wrong for PCR conditions (which again.. simulated data so this may change!!).
- Allele labels are same as ScaleHD 1.x `CAG_CAACAG_CCGCCA_CCG_CCT`.
  `42_0_1_7_2` is a loss of the CAA interruption, `19_2_1_10_2` a CAACAG duplication, etc.
  `83+_1_1_7_2` means no read spanned the entire CAG tract, so only a lower boundary call is known.
  A rough estimate of the true length is given only when the reads also limit it from
  above, and it leans on the stutter model beyond the lengths it was measured at. On
  simulated 300-base reads that happens up to about 7 CAG past the read limit (an
  estimate like 93, 89-96 for a "true" CAG90), beyond that, expect the lower boundary call alone.
  This is one element where my out-of-date science knowledge will likely result in significant
  code changes once feedback is received.
- Flags mark what deserves a manual inspection (unchanged really from previous version):

  | flag | meaning |
  |---|---|
  | `low_depth` | fewer than 500 usable molecules |
  | `low_confidence` | posterior below 0.99 |
  | `homozygous` | both alleles identical |
  | `neighbouring` | alleles one CAG apart. the hardest case to separate from stutter |
  | `close_alleles` | alleles further apart, but one's stutter is 10% or more of the other's peak, so their exact CAG are hard to tell apart. most often two long alleles |
  | `atypical` | at least one allele without the common HTT sequence structure |
  | `beyond_read_length` | an allele longer than the reads. CAG is a lower boundary call, not definitive |
  | `allele_imbalance` | one allele has under 20% of molecules |
  | `high_background` | over 5% of molecules fit neither allele |
  | `unexplained_peak` | a peak the genotype doesn't explain. third allele, contamination, mosaicism? |
  | `high_discordance` | over 15% of molecules dropped because their read base pairings disagreed |

## Checking the input went well / possible errors

The first lines `scalehd genotype` prints (and `read_outcomes` in the counts JSON)
show how reads fared:

- mostly `no_anchor` in R1 = the files are probably swapped (R1 is the CCG strand)
  or the amplicon doesn't contain the default anchors.
- many `nonconforming` = reads reach both anchors but don't fit the tract structure.
  Expect this with very poor quality or human error
- high `dropped` = read pairs often disagree. A few percent is normal sequencing error,
  much more suggests quality problems or out-of-step bases.

## Other amplicons and PCR conditions

- Other flanks. The default flanks are those of the ScaleHD 1.x reference library.
  From Python, `RepeatParser(amplicon=AmpliconSpec(...))` and `count_fastq(...,
  parser=...)` take other flanks. The CLI doesn't expose this yet i.e. TODO work
- Other stutter. Stutter priors were measured on MiSeq data from one PCR protocol
  (`scalehd.calibration.HTT_MISEQ`) with the limited data I have at the time of writing.
  Very different chemistry or cycle numbers may need a recalibrated curve, passed via
  `CallerSettings(stutter=...)`. also still WIP and subject to feedback from actual users/real data.
