"""Model-based genotype calling from per-molecule repeat structures.

A candidate genotype is a pair of alleles, possibly the same allele twice. Every
molecule in a sample is scored under a mixture::

    P(molecule) = (1 - β) · [w · M(A₁) + (1 - w) · M(A₂)] + β · background

M(A) is what PCR and sequencing make of allele A:

- the CAG length moves by stutter and somatic expansion, following the two-sided
  geometric kernel of `scalehd.calibration`, whose six ratios are fitted per
  allele around length-dependent priors;
- CCG moves by one unit either way with probability γ;
- each of the CAACAG, CCGCCA and CCT counts is misread with probability μ.

A molecule whose CAG tract ran past both reads contributes P(CAG >= its bound), and
unobserved fields are dropped out. An allele beyond read length has one unknown N >= L,
shared by all its molecules: the sample's likelihood is averaged over N in [L, L + 60),
so it must explain the data as well as a specific N would. The readable contracted
molecules then give an estimate of N. Where reads stop at L, exact candidates stay
below L, since lengths from L up can only be told apart through that average.

Near read length, many molecules' CAG tracts are read to their end, but too near the
end of the read to rule out a sequencing error etc influencing the genotype. Such an unconfirmed
count is likely to be real, but scored with appropriate caution (``unconfirmed_error``).

Only a read that ends inside the tract gives a lower genotype call boundary, and where a read
ends depends on where it starts, so lower call boundaries don't bias N.


Each candidate is fitted by maximum a posteriori. Candidates are then compared by a
Laplace approximation to their marginal likelihood, which gives a posterior probability
for the call. A homozygous candidate can only explain a real second peak with stutter
ratios its priors make implausible.
"""

from .call import call_genotype
from .candidates import candidate_alleles
from .results import SCHEMA, AlleleCall, Candidate, Flag, GenotypeCall, NoMoleculesError
from .settings import CallerSettings

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
]
