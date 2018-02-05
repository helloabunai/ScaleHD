.. _sect_genotyping:

Automated Genotyping (Gtype)
================================

The genotyping procedure of ScaleHD has two stages. The results from each stage are utilised in step with one another to act as a double-confirmation of the data. The first stage involves Digital Signal Processing (DSP). This stage's main function is to determine the structure of the Huntington Disease repeat tracts. As seen in the figure below, a typical structure consists of a CAG tract, an intervening sequence, a CCG tract and a CCT tract. In 95% of samples, the structure of the intervening sequence is one CAACAG hexamer, and one CCGCCA hexamer. This 'typical' structure is represented as: CAG_1_1_CCG_CCT. 

However, in ~5% of samples, the structure differs. These 'atypical' samples, have an intervening sequence which is different from the normally observed one. These can vary in literal size, but an example of one atypical allele would be represented as: CAG_2_1_CCG_CCT. This is important as, due to the string-similarity of the hexamers in the intervening sequence, an atypical allele can be incorrectly assigned to a typical reference with a larger CCG. For example, a sample with the genotype CAG_2_1_7_CCT would be the exact same length as CAG_1_1_9_CCT.

Thus, the DSP will scan each read in an aligned (to a typical reference) assembly, to determine the literal structure of repeat units. It does this by doing XYZ to ABC.

continue this tomorrow.....