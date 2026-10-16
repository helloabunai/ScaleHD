// ScaleHD 1.x's flags, by its InstanceReport column
// Basically just ripped straight from there
export const LEGACY_FLAGS: Record<string, { significance: string; meaning: string }> = {
  "Homozygous Haplotype": {
    significance: "Moderate/Major",
    meaning:
      "Both alleles have the same genotype. Moderate significance in individuals not affected by HD, major if genotyping HD patients or people at risk of HD.",
  },
  "Neighbouring Peaks": {
    significance: "Moderate",
    meaning: "Two alleles in the same CCG, their CAG one apart.",
  },
  "Diminished Peaks": {
    significance: "Moderate",
    meaning:
      "One peak in a CCG-homozygous sample is a large expansion with a relatively minuscule read count.",
  },
  "Novel Atypical": {
    significance: "Major",
    meaning:
      "An allele with an atypical intervening sequence, unlike the commonly observed atypical structures.",
  },
  "Alignment Warning": {
    significance: "Moderate",
    meaning: "Determining CCG returned more values than is possible (more than 2).",
  },
  "Atypical Alignment Warning": {
    significance: "Minor",
    meaning: "Re-aligning for an atypical allele, very poor quality re-alignment gave unclear data.",
  },
  "CCG Rewritten": {
    significance: "Moderate",
    meaning: "An allele's CCG (from atypical scanning) was found invalid and overwritten.",
  },
  "CCG Zygosity Rewritten": {
    significance: "Minor",
    meaning:
      "An allele deemed CCG-heterozygous (typical reference) was found to be an atypical CCG-homozygous allele.",
  },
  "CCG Uncertainty": {
    significance: "",
    meaning:
      "From the legacy code: its CCG peak checks were ambiguous, e.g. forward and reverse reads showed different numbers of CCG peaks.",
  },
  "CCT Uncertainty": {
    significance: "Minor",
    meaning:
      "Atypical scanning found no clear agreement on the CCT length (e.g. CCT2 55%, CCT3 45%).",
  },
  "SVM Failure": {
    significance: "Major",
    meaning: "The SVM's confusion matrix was inconclusive, and CCG zygosity had to be bootstrapped.",
  },
  "Differential Confusion": {
    significance: "Major",
    meaning:
      "The allele sorting was confused between a potential neighbouring peak or a homozygous haplotype.",
  },
  "Missed Expansion": {
    significance: "Major",
    meaning:
      "One genotyping algorithm believes an allele expanded with low reads, the other believes otherwise.",
  },
  "Heuristic Filtering Success": {
    significance: "",
    meaning:
      "From the legacy code: atypical scanning's heuristics chose the second allele. A pass, not a warning.",
  },
  "Peak Inspection Warning": {
    significance: "Minor",
    meaning:
      'At least one allele failed the minimum read count inspection. Common in "bad" sequencing data.',
  },
  "Low Distribution Reads": {
    significance: "Dependent",
    meaning:
      "At least one allele's CAG read distribution (200 forward reads long) has notably few reads.",
  },
  "Low Peak Reads": {
    significance: "Major",
    meaning: "An allele's peak has very few reads (200 or fewer). Genotyping is hard here.",
  },
};

// The stage 1.x failed in, by its own name for it (its InstanceReport's Exception Raised).
export const LEGACY_FAILED_STAGES: Record<string, string> = {
  SeqQC: "quality control",
  SeqALN: "alignment",
  DSP: "atypical scanning",
  "SeqRE-ALN": "realignment",
  Genotype: "genotyping",
  SNPCalling: "SNP calling",
  "Report/Graph": "writing its report",
};

/** 1.x's guide to its confidence (Definitions.rst). */
export function legacyConfidenceInWords(confidence: number): string {
  if (confidence > 80) return "good";
  if (confidence > 60) return "believable";
  if (confidence < 55) return "of major significance: may justify manual inspection";
  return "may justify manual inspection";
}
