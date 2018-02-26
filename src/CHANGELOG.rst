.. _sect_changelog:

Version 0.300
-------------

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