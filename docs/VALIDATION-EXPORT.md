# Validation export

Table of numbers and magic from the model that was used during genotyping.
Shareable (maybe? pending confirmation) with developers to improve model after tests on real data,
as raw data files may be unavailable to devs.
Download it from a job's page ("Validation export (CSV)") once the job has
finished or been cancelled, or from `GET /api/jobs/{id}/validation.csv`.

## What's in it

A row per sample:

- Validation ID and batch. The ID is the job's number and the sample's place in it
  (`7-003`). The batch is the folder the reads came from, as `B1`, `B2` and so on in order.
- Columns for the lab to fill in (left blank): lab genotype, the method it came from,
  replicate group, sample type and notes. See below. All may not be required.
- The call. status, genotype, quality, posterior and flags, plus the three runner-up
  genotypes and their posteriors.
- The fitted values. Background, CCG slippage and misread rate, any unexplained peaks,
  and for each allele its structure, CAG (with a `+` and a rough estimate when it's too long
  to read), CCG, share and molecules, its six fitted stutter ratios, and the ScaleHD 1.x style
  backward slippage, somatic mosaicism and expansion and contraction indices.
- Each allele's molecules around its peak. Molecules of the allele's structure at each
  CAG length from N-20 to N+30 (the reach of the stutter model), space separated. Both
  alleles' molecules are in there when they share a structure. Blank for an allele too long
  to read.
- Read quality. Molecules, complete, partial, dropped and unusable, each read's
  outcomes (spanning, truncated, no anchor, nonconforming), which tracts the two reads
  disagreed on, and the read length.
- Truth and whether it matched. For simulated samples only.
- The kind of error (if present), for a failed sample (e.g. `FileNotFoundError`), without its specific message.
- The ScaleHD and server versions that ran the job, the method and every other setting.

## What's left out

- Sample names, file paths and the job's name.
- Error messages, which can name a sample or a file, only the kind of error is kept.

The genotypes are still there. A Huntington's disease genotype is sensitive health
information, and in a small cohort a genotype with a few other columns could point to a
person. Whether this table can be shared is a question for whoever handles data governance
where the data is held so let's see if evne this can be shared

## Matching rows back to samples

The ordinary results export ("Export results (CSV)") has the sample names and the same
Validation ID column. Keep that one in the lab and share only the validation export.

## Filling in the lab columns

- Lab genotype: the genotype from another method, in either form, `17_1_1_7_2/43_1_1_7_2`
  or just the CAG lengths `17/43`.
- Lab genotype method: optional i guess but e.g. capillary sizing, ScaleHD 1.x, manual review, Sanger.
- Replicate group: also optional but the same label on repeats of one sample, so they can be compared.
- Sample type: optional, blood, brain, cell line and so on (somatic expansion differs by tissue).
- Lab notes: optional, if any, anything else worth knowing about the sample or the run.
