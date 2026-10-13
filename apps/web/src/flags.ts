// What each genotype flag means, as in the README's flag table.
export const FLAGS: Record<string, string> = {
  low_depth: "Fewer usable molecules than the minimum (500 by default).",
  low_confidence: "The posterior is below the threshold (0.99 by default).",
  homozygous: "Both alleles are identical.",
  neighbouring: "The alleles are one CAG apart. The hardest case to tell from stutter.",
  close_alleles:
    "The alleles are close enough that one's stutter is a fair share (10% or more) of the other's peak, so their exact CAG are hard to tell apart. Most often two long alleles.",
  atypical: "At least one allele without the common HTT sequence structure.",
  beyond_read_length: "An allele longer than the reads. Stated CAG is presumed as lower bound.",
  allele_imbalance: "One called allele has under 20% of the molecules.",
  high_background: "Over 5% of molecules fit neither allele.",
  unexplained_peak:
    "A peak the genotype doesn't explain. A third allele, contamination, mosaicism?",
  high_discordance: "Over 15% of molecules were dropped because their read base pairings disagreed.",
};
