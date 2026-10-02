// What each genotype flag means, as in the README's flag table.
export const FLAGS: Record<string, string> = {
  low_depth: "Fewer usable molecules than the minimum (500 by default).",
  low_confidence: "The posterior is below the threshold (0.99 by default).",
  homozygous: "Both alleles are identical.",
  neighbouring: "The alleles are one CAG apart. The hardest case to tell from stutter.",
  atypical: "An allele without the common 1_1_x_2 intervening and CCT structure.",
  beyond_read_length: "An allele longer than the reads. Stated CAG is presumed as lower bound.",
  allele_imbalance: "One called allele has under 20% of the molecules.",
  high_background: "Over 5% of molecules fit neither allele.",
  unexplained_peak:
    "A peak the genotype doesn't explain. A third allele, contamination, mosaicism?",
  high_discordance: "Over 15% of molecules were dropped because their read base pairings disagreed.",
};
