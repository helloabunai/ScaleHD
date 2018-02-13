.. _sect_definitions:

Miscellanious Definitions
================================

Within the instance-wide output table produced by ScaleHD, there are many flags or data entities which require explanation. Throughout the development of ScaleHD, we ended up determining a range of characteristics that indicate the believability of data produced from amplicon sequencing when attempting genotyping via our multiple reference library-based method. These characteristics provide us with a heuristic method to determine if an attempt at automated genotyping was successful, or not. Here, we define these characteristics and what they mean in literal terms. They range in importance, but all are useful in creating a representation of data quality. Some are self-explanatory, but are explained anyway.

The significance levels are as follows:

 * N/A -- This means the flag contains discrete information and does not need to be interpreted in regards to genotyping quality.
 * Dependent -- This flag may be significant, depending on other flags. For example, a high level of somatic mosaicism may be an indicator of poor genotyping quality when the CAG repeat tract size is within the non HD-causing allele size range.
 * Minor -- This entity is of minor significance and in the vast majority of samples will not be a deterministic factor for genotyping quality.
 * Moderate -- This entity is of moderate significance. It is unlikely to render a sample's genotype invalid on its own but may contribute to inaccurate genotyping.
 * Major -- This entity is of major significance and is strongly associated with genotyping quality. If any major informative flags are raised, it is recommended to manually inspect this sample.

+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| ScaleHD Flag               | Significance | Definition                                                                                                             |
+============================+==============+========================================================================================================================+
| SampleName                 | N/A          | Literally the sample name/label taken from the file                                                                    |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| Pri/Sec GTYPE              | N/A          | The genotype returned for each allele (e.g. CAG_1_1_CCG_2)                                                             |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| Status                     | Dependant    | Whether the structure of the intervening sequence is typical (1_1) or atypical (anything else).                        |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| BSlippage                  | Dependent    | The amount of backwards slippage, relative to each allele's peak. Calculated as [(n-1 to n-5) / n].*                   |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| Somatic Mosaicism          | Dependent    | The amount of somatic mosaicism, relative to each allele's peak. Calculated as [(n+1 to n+10) / n].*                   |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| Confidence                 | Major        | The percentage confidence in a genotype call. See the confidence calculation subsection for info.                      |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| Exception Raised           | N/A          | If the pipeline reached a fatal error on a sample, the stage at which it crashed is listed here.                       |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| Homozygous Haplotype       | Moderate     | If a sample has a homozygous haplotype -- i.e. both alleles have the same genotype.                                    |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| Neighbouring Peaks         | Moderate     | Two alleles exist within the same CCG dimension, with CAG values being separated by 1 value.                           |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| Diminished Peaks           | Moderate     | One peak in a CCG-homozygous sample is a large expansion with a relatively minuscule read count.                       |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| Novel Atypical             | Major        | An allele with an atypical intervening sequence, different to that of the commonly observed atypical structures.       |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| Alignment Warning          | Moderate     | When determining CCG values, more values were returned than is possible (i.e. more than 2 results).                    |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| Atypical Alignment Warning | Minor        | When re-aligning for an atypical allele, particularly awful quality re-alignment produced unclear data.                |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| CCG Rewritten              | Moderate     | An allele's CCG value (from DSP) was determined invalid and overwritten.                                               |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| CCG Zygosity Rewritten     | Minor        | An allele was (typical reference) deemed CCG-heterozygous, but detected to be an atypical CCG-homozygous allele.       |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| CCT Uncertainty            | Minor        | When using DSP to determine the CCT tract length, there was no clear agreement (e.g. CCT2 = 55%, CCT3 = 45%).          |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| SVM Failure                | Major        | The confusion matrix produced by the SVM was inconclusive, and CCG zygosity had to be bootstrapped.                    |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| Differential Confusion     | Major        | The allele sorting algorithm was confused between a potential neighbouring peak, or homozygous haplotype.              |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| Peak Inspection Warning    | Minor        | At least one allele failed minimum read-count distribution threshold inspection. Common in "bad" sequencing data.      |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| Low Distribution Reads     | Dependent    | At least one allele's CAG read distribution (FW-200 length) contains a noteworthly low number of reads.                |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
| Low Peak Reads             | Major        | In a given allele's read distribution, the n value contains a very low number** of reads. Genotyping is hard, here.    |
+----------------------------+--------------+------------------------------------------------------------------------------------------------------------------------+
``
*n denotes the number of reads for the modal allele
**very low reads is defined as an n value containing <=200 reads
``

Confidence Calculation
~~~~~~~~~~~~~~~~~~~~~~

For each allele, ScaleHD calculates the confidence level in the provided genotyping result. This information is taken from a variety of sources, and attempts to paint an evidence-based picture of the data quality, and resultant genotype confidence. Each allele starts with 100% confidence, and penalties are applied when certain data characteristics were discovered throughout the genotyping process. Follows is a list of evidence used to best determine each allele's confidence level:

 * If the First Order Differential peak confirmation stage required to re-run itself, with a lower threshold. More re-calls results in a higher penalty.
 * Rare characteristics, such as homozygous haplotypes, or neighbouring/diminished peaks, incur a penalty.
 * Atypical alleles are treated with more caution, and scores are weighted slightly more severely than typical alleles.
 * Simple data aspects such as total read count within a sample/distribution/peak are used.
 * Mapping percentages are taken into account, albeit as a minor factor within this algorithm.
 * "Fatal" errors, such as Differential Confusion, incur a significant penalty.

Any confidence score is capped at 100%. If the quality of data in a particular sample is high enough for alleles to be awarded a confidence score higher than 100%, they are reported as 100%, regardless. Generally, a 'good' score is anything over 80%, and we have found that samples returning a score of over 60% are considered believable. Anything less than this may justify manual inspection.