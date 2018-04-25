.. _sect_changelog:

Version 0.310
-------------

This is a minor update to ScaleHD. SNP calling implementation has begun, but requires polishing,.

 * Fixed a bug where genotyping would complete, but raise an exception at the end of the genotyping module, due to particular arrays not being flattened.
 * Implemented Picard/GATK into the SNP calling module of ScaleHD.
 * Updated the requirements for Picard/GATK to be integrated with ScaleHD on the user's system.
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
