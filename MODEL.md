# The new genotyping model

This section explains the "New (model-based)" genotyping method i.e. how it turns a sample's
molecule count into two called alleles, a confidence in the call, and analysis flags. The code is in
[`packages/core/src/scalehd/genotype/model/`](packages/core/src/scalehd/genotype/model/),
and every default below is in its `settings.py`. Users will eventually be able to fine-tune parameters
from the web interface, but this is not implemented yet, because the model may already undergo
significant changes when access to real data is received.

The PCR slippage/stutter assumptions come from training data used by ScaleHD 1.x for a support vector machine,
and other algorithm components have only been tested on simulated reads so far.
Expect the numbers, and maybe some of the model, to change once real data is available.

## What the inputs are

As we are not using sequence alignment for this process, we instead use raw read counts from the sequencing
data.

Counting (`scalehd count`) reads every molecule's repeat structure, its five counts
`CAG_CAACAG_CCGCCA_CCG_CCT`, straight from its reads. Each count has a status:

| status | meaning |
|---|---|
| exact | the tract was read from end to end |
| lower bound | the read ended inside the tract, so it is at least this long (e.g. `83+`) |
| unconfirmed | the tract's end was seen, but too near the end of the read to rule out e.g. sequencing error (`80~`) |
| unobserved | the repeat tract was not seen in the read|

The caller works with these distinct observations and how many molecules had each.

## Candidate alleles and genotypes

The genotype caller attempts fitting candidate alleles via:

- the 6 most common structures among molecules read in full (counting a CAG that is
  unconfirmed but otherwise complete as read in full),
- CAG +/- 1 of the top 3, as PCR slippage/stutter can make the true allele the second most common
  peak (backwards and/or forward slippage),
- one allele beyond read length, when enough molecules ran past both reads (at least
  20, and at least 3% of the sample). Its CAG is `L`, the most common lower bound, and its
  other counts are the most common exact values among those molecules. Exact candidates
  at `L` or above are dropped. From `L` and up, lengths can only be told apart through the
  averaging described below.

A candidate genotype is any pair of these scenarios, including the possibility of the same allele twice
(homozygous).

## The likelihood of one molecule

Every molecule `x` is scored under a mixture of the genotype's alleles `A1` and `A2`, plus a
background:

```
P(x) = (1 - β) · [ w · M(x | A1) + (1 - w) · M(x | A2) ] + β · B(x)
```

| symbol | meaning |
|---|---|
| `w` | allele balance: the shorter allele's share of molecules (1 when homozygous) |
| `β` | background: the share of molecules that come from neither allele |
| `B(x)` | background spread evenly over every structure: 250 CAG values, 4 CAACAG, 4 CCGCCA, 25 CCG and 5 CCT |
| `M(x | A)` | what PCR and sequencing make of allele `A` |

`M` multiplies one term per count. A count that wasn't observed contributes nothing:

```
M(x | A) = K(x_CAG | A) · C(x_CCG | A) · T(x_CAACAG | A) · T(x_CCGCCA | A) · T(x_CCT | A)
```

### CAG and the slippage/stutter kernel

Again, this was derived from limited access to real data, and is likely to change.
PCR slippage spreads an allele of `N` CAG over nearby lengths. The heights of the peaks
relative to `N` follow six ratios, three on each side:

```
h(0)  = 1
h(-1) = c1                         c1 = (N-1)/N
h(-k) = c1 · c2 · ct^(k-2)         c2 = (N-2)/(N-1), ct = each further step down
h(+1) = e1                         e1 = (N+1)/N
h(+k) = e1 · e2 · et^(k-2)         e2 = (N+2)/(N+1), et = each further step up
```

The first two steps on each side get their own ratios because real expanded alleles have
a sharp `N+1` step followed by a long, flatter somatic tail, which one geometric decay
cannot represent. In blood DNA, the expansion side includes ordinary somatic expansion.

Heights are zero beyond 20 below `N` and 30 above, and below CAG 1. Normalised, they give
the kernel:

```
K(x | N) = h(x - N) / Z(N)        Z(N) = the sum of h over every length the kernel reaches
```

A small floor `φ` spreads a share of the allele's molecules evenly over CAG 1 to `S`,
where `S` is the longest CAG seen (at least 40) plus 10. Without it, a long flat tail
would distort the stutter fit and shift `N`:

```
K'(x | N) = (1 - φ) · K(x | N) + φ / S
```

A CAG that is only a lower boundary call (e.g. `83+`) `b` scores the chance of a tract at least that long,
`P(CAG >= b | N)`: the sum of `K'` from `b` up.

An unconfirmed count `x` is probably right, but the tract may go on past a sequencing
error. With `ε = 0.03`:

```
(1 - ε) · K'(x | N) + ε · P(CAG >= x + 1 | N)
```

### CCG

```
C(x | A) = 1 - 2γ - ψ + ψ/25     x is the allele's CCG
           γ + ψ/25              x is one unit away
           ψ/25                  x is two or more units away
           γ^d                   x is a lower bound d units above the allele's CCG
           1                     x is a lower bound at or below the allele's CCG
```

`γ` is CCG slippage by one unit each way. `ψ` is a CCG read as any value at all: in
ScaleHD 1.x alignments, a few percent of reads that didn't span the CCG tract landed on
the wrong CCG entirely, and without this term they distort the other allele's fit.

CCG was typically a lot cleaner than CAG in sequencing, from memory, but assumptions will
be naturally corrected if required.

### CAACAG, CCGCCA and CCT

Each of these counts is misread with probability `μ`:

```
T(x | A) = 1 - μ      x is the allele's count
           μ / 3      x is a different count, or a lower bound above the allele's
           1          x is a lower bound at or below the allele's
```

## The sample's likelihood

The sample's log-likelihood sums every distinct observation's `log P(x)`, times the
number of molecules that had that observation. Reads are PCR copies of a limited number of input
templates, not independent molecules, so the counts are scaled down to at most 3,000
molecules in total. Without that, at 10^5 reads tiny misfits in peak shape would outweigh
every prior assumption.

### An allele beyond read length

All the molecules of an allele beyond read length come from one template length `N`,
which is unknown but at least `L`. The model doesn't fit `N`. It averages the whole
sample's likelihood over the 60 lengths from `L` to `L + 59`:

```
likelihood = (1/60) · Σ over N in L .. L+59 of  Π over molecules x of  P(x | N)
```

So the allele has to explain the data across that range, rather than at whichever single
`N` happens to fit best. In comparison to ScaleHD 1.x's logic, this approach of asking the
explanation of a genotype observation to fit both in local and "global" terms, provides
(hopefully) an automated genotyping tool which is a lot more robust and performant.

## Parameters and priors

Each genotype has these parameters:

- six stutter ratios per allele,
- the balance `w` (heterozygous only),
- the shared `γ`, `ψ`, `μ`, `φ` and `β`.

All are fitted on log or logit scales, each with a Gaussian prior. The priors are given
as their medians:

| parameter | prior median | prior SD (log or logit scale) |
|---|---|---|
| stutter ratios | the `HTT_MISEQ` curve at the allele's CAG | 0.35, 0.45, 0.8, 0.6, 1.0, 1.2 |
| `w` balance | 0.5 | 0.6 |
| `γ` CCG slippage | 0.01 | 1.0 |
| `ψ` CCG read as anything | 0.001 | 1.5 |
| `μ` CAACAG/CCGCCA/CCT misread | 0.002 | 1.0 |
| `φ` floor | 0.005 | 1.5 |
| `β` background | 0.002 | 1.5 |

The `HTT_MISEQ` curve (`scalehd/calibration.py`) is made from per-length medians of the
six ratios from training data used by ScaleHD 1.x. That is ~600 MiSeq samples, leaving out
alleles within 8 CAG of a partner with the same CCG, and alleles with fewer than 500
reads at `N`. The curve is interpolated between the measured lengths and held constant
beyond them. Longer alleles stutter more, so the assumption moves with the allele.

This is where the model may be incomplete, and more data will allow me to improve this.

`c1` and `e1` are capped at 1, so an allele is always the tallest peak of its own
molecules, as in the usual sizing convention.

## Fitting and comparing genotypes

Fitting a genotype means finding the parameters that maximise likelihood times prior
(L-BFGS-B, within bounds). There are three stages:

1. Screen: Every candidate genotype is scored with low-effort. Its parameters stay at their
   prior medians, apart from the balance, which is optimised. Each score is penalised
   `0.5 · k · log(n)` for its `k` parameters and `n` molecules (after the 3,000 cap).
2. Fit: The 6 best are fitted in full, together with the best homozygous and the
   best heterozygous genotype if they aren't among the 6.
3. Confirm: All other candidates are fitted again, starting from the leader's background, floor
   and misread values. Those parameters can settle in different places for different
   candidates, which would leave the comparison down to luck.

Genotypes are compared by their marginal likelihood (the evidence for them, averaged
over their parameters), using a Laplace approximation around the fit `θ̂`:

```
score = log L(θ̂) + log prior(θ̂) - 0.5 · log det H(θ̂)       (constants dropped)
```

`H` is the curvature of minus the log posterior at the fit. Its eigen values are floored at
the widest prior's curvature, since the posterior can't be flatter than the prior. Each
genotype's posterior probability is the share of `exp(score)` over all the fitted
genotypes.

This is what stops a homozygous genotype from explaining a real second peak. It can
only do so with stutter ratios far from their priors, which this approach penalises.

## Settling each allele's exact CAG

The full fit decides which alleles are there. Each allele's exact `N` is then decided
from its own peak alone, so that distant shoulders and junk can't move it:

- The molecules used are those with the allele's structure and a CAG within +/- 6 of `N`.
  The other allele, the floor and the background stay as fitted.
- Each `N` from 2 below to 2 above is tried as an exact allele, with its stutter refitted
  and Laplace-scored. Stutter trades off against `N`, so holding it fixed would claim `N`
  way too precisely.
- If over 5% of the peak region wasn't read in full (near read length), that window cuts
  across the molecules that place `N`. Every molecule is used instead.
- `N` stays below `L` when the sample has an allele beyond read length.
- This step is skipped for an allele beyond read length, and for alleles within 3 CAG of
  another allele with the same structure. Their peaks overlap, so the joint fit decides
  them.

If the best `N` changed, the genotype is fitted again with it.

## Confidence and quality

```
P(call) = P(this pair of alleles) · Π over alleles of P(its exact N)
```

`P(this pair of alleles)` adds up every fitted genotype with the same alleles apart from
CAG differences of at most 2. The other `N` values tried, and the other genotypes, are
the alternatives. The top 3 are reported.

```
quality = -10 · log10(error)        error = the larger of 1 - P(call) and the alternatives' total
```

The quality is capped at 99 (nobody's perfect 💅🏻), so e.g. 20 means 1 in 100 wrong. It covers stutter 
and sampling noise. It doesn't cover the model itself being wrong for a sample's PCR conditions (yet?).

## An estimate for an allele beyond read length

For an allele beyond read length, the caller tries each `N` from `L` to `L + 59` as an
exact allele, refitting its stutter each time. It tries every 5th `N` first, then every
`N` where those leave any chance.

If the chance at `L + 59` is still above 1% of the peak's, the reads set no upper limit
and no estimate is given. Any range would only show how far the search went. Otherwise
the estimate is the most likely `N`, with the 5% to 95% range, e.g. 93 (89-96) for a true
90 on simulated 300-base reads.

## Per-allele figures and flags

Each molecule is shared out between the alleles by how likely each made it. That gives
each allele its molecules and its own CAG distribution. From those come the ScaleHD 1.x
style ratios:

- backward slippage: molecules at `N-2` and `N-1`, over those at `N`,
- somatic mosaicism: `N+1` to `N+10`, over `N`,
- the expansion and contraction indices: the mean shift above and below `N`.

A peak of molecules read in full, holding at least 2% of the sample and well beyond
what the fit expects there, is reported as unexplained. The flags, and their thresholds,
are listed in [USING-FASTQ.md](USING-FASTQ.md#reading-the-result).

## How well it performs (estimates)

Most testing of this model is partly circular due to lack of data (at the moment). The simulator stutters
with the same curve the model's priors come from, and `legacy_matrix.py` scores the caller on the matrix
those priors were measured on. Three benchmarks in `packages/core/benchmarks/` test what
flaws are maybe hidden by lack of data. The results below are from 2026-10-05, on simulated 300-base reads.

### When stutter isn't what the caller expects

`stutter_robustness.py` simulates the 19 scenarios of `genotype_simulated.py`, 3 seeds
each, with stutter different to the model priors, and calls them with the priors unchanged.

| simulated stutter | right | wrong at Q 20 or more | wrong calls the posteriors expected |
|---|---|---|---|
| as calibrated | 57/57 | 0 | 1.4 |
| half the stutter | 57/57 | 0 | 0.2 |
| double the stutter (N-1/N capped at 0.95) | 56/57 | 0 | 2.7 |
| longer tails | 48/57 | 1 | 4.2 |
| shorter tails | 57/57 | 0 | 0.0 |
| as if 10 CAG longer | 56/57 | 0 | 1.4 |
| as if 10 CAG shorter | 56/57 | 0 | 1.9 |
| drawn per sample, within the prior spread | 53/57 | 0 | 3.0 |
| drawn per sample, within twice the prior spread | 45/57 | 4 | 1.5 |

- Inside the priors' range, wrong calls are rare, and each came flagged `low_confidence` (good).
- With longer tails, or stutter at twice the prior spread, long alleles (CAG 55 to 80)
  are called one CAG off. Most are flagged `low_confidence`. Four weren't: 80 called 79
  (Q 20), 55 called 56 (Q 30), and 70 called 71 twice (Q 31 and 38) (not great not terrible).
- A homozygous 21 was called 20/21 at Q 32, but flagged `neighbouring` (which is fine i guess).

### Whether the confidence is "honest"

`posterior_calibration.py` simulates 200 random samples three times. Genotypes via
`scalehd simulate-run` draws them, plus neighbouring and near alleles, at 150 to 5,000
read pairs. If the posteriors are "honest", a call at 0.99 is wrong about 1 time in 100.

| simulated stutter | wrong | wrong calls the posteriors expected | right, among calls at 0.999 or more |
|---|---|---|---|
| as calibrated | 3 | 4.8 | 166/166 |
| drawn per sample, within the prior spread | 13 | 6.8 | 151/155 |
| drawn per sample, within twice the prior spread | 29 | 10.3 | 150/153 |

- With stutter as calibrated, the posteriors are "honest", if slightly cautious (arguably good).
- When each sample's stutter varies as much as the priors themselves allow, the
  posteriors are about twice as confident as they should be, and some calls at 0.999 or
  more are wrong.
- The wrong ones were almost all two alleles a few CAG apart, both 37 or more:
  - 44/45 called 44/46, and 43/45 called 44/46,
  - at twice the spread, also 40/41 called 41/49, and 37/38 called 38/41,
  - plus one long allele, 57 called 58.
- None of these were flagged (bad). 37/38 called 38/41 would move an allele from the
  reduced-penetrance range into full penetrance (clinical significance).

### Priors measured on the matrix they're scored on

`prior_crossval.py` rebuilds the stutter curve from a random half of the training matrix.
It then genotypes every sample twice, with the curve from the half it isn't in, and with
the curve from its own half.

- Exact calls: 407/594 (68.5%) with the curve that never saw the sample, 408/594 (68.7%)
  with the one that did. One sample changed.
- So the legacy matrix score isn't flattered by the priors coming from the same matrix.

### What this means

- When its stutter assumptions hold, the model calls well and its confidence can be
  trusted.
- How far real PCR strays from those assumptions is unknown until there's real data. If
  real samples vary as much as the priors allow, the quality is optimistic, mostly for:
  - the exact CAG of long alleles (one off),
  - two alleles a few CAG apart that are both expanded or close to it.
- Ideas, none tried yet:
  - a flag for two alleles within about 3 CAG of each other when both are above about 35
    (`neighbouring` only covers 1 apart),
  - wider, or heavier-tailed, stutter priors for long alleles,
  - letting the exact-`N` step carry more of the trade-off between stutter and `N`,
  - real samples with known genotypes, to set the prior spread from data rather than
    from ScaleHD 1.x alignments.

To rerun them:

```sh
uv run python packages/core/benchmarks/stutter_robustness.py      # about 9 minutes on AMD 5950x cpu
uv run python packages/core/benchmarks/posterior_calibration.py   # about 6 minutes on AMD 5950x cpu
uv run python packages/core/benchmarks/prior_crossval.py          # about 16 minutes on AMD 5950x cpu
```
