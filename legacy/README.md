ScaleHD: Automated Huntington Disease genotyping with end-to-end allele masking
===============================================================================
ScaleHD is a package for automating the process of genotyping microsatellite repeats in Huntington Disease data.
We utilise machine learning approaches to take into account natural data 'artefacts', such as PCR slippage and somatic
mosaicism, when processing data. This provides the end-user with a simple to use platform which can robustly predict genotypes from input data.

By default, input is a pair of unaligned .fastq sequence data -- both forward and reverse reads, per sample. We utilise both forward and reverse
reads in order to reduce the complex dimensionality issue posed by Huntington Disease's multiple repeat tract genetic structure. Reverse reads allow
us to determine the current sample's CCG state -- this provides us with a mechanism by which to more easily call the entire genotype. Forward reads
are utilised in a similar approach, to determine the CAG and intervening structure.

The general overview of the application is as follows:
1) Input FastQ files are subsampled, if an overwhelming number of reads are present. This can be overruled with the -b flag.
2) Sequence quality control is carried out per the user's instructions. We reccomend trimming of any 5-prime spacer+primer combinations, for optimal alignment.
3) Alignment of these files, to a typical HD structure (CAG_1_1_CCG_2) reference, is carried out.
4) Assemblies are scanned with Digital Signal Processing to detect any possible atypical structures (e.g. CAG_2_1_CCG_3).
4.1) If no atypical alleles are detected, proceed as normal.
4.2) If atypical alleles are detected, a custom reference is generated, and re-alignment to this is carried out.
5) With the appropriate allele information and sequence assembly(ies) present, sampled are genotyped.
6) Output is written for the current sample; the procedure is repeated for the next sample in the queue (if present).


Check the full documentation at http://scalehd.rtfd.io
