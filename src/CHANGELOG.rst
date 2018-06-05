.. _sect_changelog:

Verson 0.316
------------

* Added some minor documentation for SNP Calling (_sect_genotyping)
* Heuristic allele filtering engine has been completely rewritten to not be absolute garbage.
* Parallelised the DSP module within ScaleHD to execute on multiple contigs of data at once, if enabled.
* Parallelisation introduced issue with allele structure incrementing objects would behave improperly -- this is now fixed.
* Disabled subsampling of aligned assemblies (due to multi-threading speedup; no longer required).
* Implemented broad error catching around SNP calling libraries, instead of just exiting upon failure.
* Fixed bug with PDF rendering of result distributions utilising an incorrect value for aligned read counts.
* Fixed bug where atypical alleles which changed from CCG-homozygous to CCG-heterozygous were not identified.
* Fixed error where the heuristic filtering engine suspects an expanded allele, but ended up calling a homozygous haplotype.
* Casting issue where two alleles returned different dimension-shaped arrays for FOD genotype calling, was resolved.

Version 0.314/5
---------------

 * Fixed homozygous haplotype casting error
 * Fixed diminished alleles being skipped (or not flagged) in particular cases of read drop-off in homozygous expansions

Version 0.313
-------------

 * Fixed a rare error wherein graphs would not be rendered where an atypical allele rewrote the CCG-zygosity from heterozygous to homozygous.
 * Added a flag for when the two core genotyping algorithms cannot agree on the status of one allele; this manifests as an expanded allele being missed due to significantly low read count.
 * Allele sorting algorithm has been tweaked to correct some mistakes in my garbage code.
 * Fixed rare error where FastQC would be executed on incorrect data.
 * Fixed certain genotyping flags being applied on a sample wide basis as opposed to an individual allele basis.

Version 0.312
-------------

 * Added an additional (optional) pre-processing stage, including sequence demultiplexing via Batchadapt.
 * CCG First order differential bugfix in situations where peak-calling returned multiple variables when unexpected.
 * Added Batchadapt to the required python package list for ScaleHD. Installed automatically from PIP where possible.

Version 0.311
-------------

 * Moron hotfix for dumb reverse aggregate distribution bug I introduced with v0.310

Version 0.310
-------------

This is a minor update to ScaleHD. SNP calling implementation is now in alpha.

 * Fixed a bug where genotyping would complete, but raise an exception at the end of the genotyping module, due to particular arrays not being flattened.
 * Implemented Picard/GATK/Freebayes into the SNP calling module of ScaleHD.
 * Added PyVCF as a Python library requirement for scraping data from variant calls.
 * Modified the requirements for Picard/GATK to be integrated with ScaleHD on the user's system $PATH.
 * Added Freebayes to the list of required binaries in __backend; addition user $PATH check
 * Added new XML flag for user to specify a strictness value, for determining legitimate SNP calls.
 * Minor codebase re-arranging in preparation for Digital Signal Processing to be replaced by a c++ binary, for performance.

Version 0.300
-------------

We now consider version 0.300 a "release-candidate alpha", if such a thing exists. I.E. The functionality performs as desired, 99% of the time (figure not accurate and i am not legally liable for any repercussions of assuming ScaleHD is 99% accurate haHAa). From this point onwards, new releases will contain new features, or a large collection of bug fixes. Minor iterations are (hopefully) over.

 * Removed Rpy2 and R-interface codebase in preparation for switching bayesian confirmation model to a native python library.
 * Added additional flag for ScaleHD output, describing how many reads that mapped to multiple references were removed (if enabled by the user).
 * Switched output rendering pipeline from Prettyplotlib to Seaborn (PPL is no longer supported).
 * Minor backend modifications in relation to the above.
 * SKLearn deprecation on label encoder fixes
 * Minor genotyping fixes (thresholds)

Version 0.252
-------------

 * Modified the N-Aligned distribution logic to utilise pre-smoothing data distribution as opposed to post-smoothing.
 * Bugfix with label in (a)typical allele being assigned an estimated CAG attribute which was not an integer.
 * FastQ subsampling workflow modified to remove possibility of incorrect percentages applying to genotyping confidence.
 * Fixed the algorithm which calculates Somatic Mosaicism for each allele (i.e. no longer reading from incorrect attributes).
 * Some other stuff that I forgot.

Version 0.251
-------------

 * Removed the redundant workflow codebase for Assembly processing (i.e. using BAM as input; feature not required/desired anymore).
 * Refactored the input method that the user can specify to subsample input reads, or not.
 * Scope fix for instances that do not use SeqQC.
 * Alternative shell pathing check for requisite binaries fix (e.g. using zsh instead of bash)

Version 0.250
-------------

 * CCG distribution cleanup threshold tweaks
 * Added handler for atypical-typical 50:50 read ratio assembly contigs.
 * Added a threshold context manager for Neighbouring Allele Peak algorithm.
 * Added differential confusion flag for samples which ScaleHD cannot sort via heuristics.
 * Begun to implement Polymorphism detection..
