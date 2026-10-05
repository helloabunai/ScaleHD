# What's different to ScaleHD 1.0x

The original implementation of ScaleHD was a command line python package, which took settings
from users via an XML file, and executed pipeline steps based on what was requested. It was a 
linear pipeline that didn't tolerate unknown situations very well. Many aspects could have been
designed in a more robust way but everyone starts somewhere so it is what it is.

It would align sequencing reads with BWA-MEM to 4,000 synthetic references,
scanned the presented structure for atypical alleles, and re-aligned to custom references when
any were found. Genotyping was done with fuzzy logic that in hindsight could be improved
a lot.

This re-write of the entire software (interfaces, algorithms, everything) to be a docker container 
(for ease of shipping), which provides a web based front-end, where users can specify settings,
view past jobs/results, run specific ScaleHD features based on requirements, and export results 
to PDF files (or similar). While there will be an API for calling backend functionality from the
web interface, the backend/python package can also be used as before, i.e. a command line
interface, if users prefer.

By default ScaleHD now uses a statistical model based approach for genotyping data without doing
any sequence alignment. This results in massively faster genotyping times because the majority of 
ScaleHD 1.x analysis time was spent doing sequence alignment. This model is a first-implementation
approach and may change massively - data used for testing has been generated with a simulator that
I have also written to generate HTT FastQ files based on my memory and results in papers that
the Monckton group published during my time working with them. This likely means the first implementation
will get things wrong when presented with real data, but the structure exists so that modifications
to the statistical model in the future will be straight forward.

Comparisons between ScaleHD 1.x and this version are thus a concern for scientific reproducability and as
a result I will eventually incorporate the exact same code/algorithm that ScaleHD 1.x uses for genotyping.
Users will be able to choose which genotyping approach they want to use when configuring a job to be run.

This re-write, or ScaleHD 2, or whatever this is called, reads the repeat structure straight
from each read instead:

1. Short sequences at the flank ends next to the repeat locate it in the
   read. Primers, spacers, adapters and off-target reads should thus need no trimming first.
2. The bases between the anchors are split into the established Huntington structure
   `(CAG)n (CAACAG)a (CCGCCA)b (CCG)m (CCT)k`, tolerating substitutions and single-base
   insertions or deletions.
3. R1 and R2 are the same molecule. When both see a repeat tract and
   disagree, the molecule is dropped, because sequencing errors rarely coincide while
   PCR stutter is shared. For long alleles, R1 supplies the CAG tract and R2 the CCG
   and CCT tracts.
4. A read that ends inside the repeat only gives a lower boundary call (e.g. `83+`).
   It is never forced onto a reference with confidence. A read that ends just after the end of a repeat tract
   gives that tract's count, but marked as an estimate (`80~`). A sequencing error near the end of
   a read could create misleading data (even if this is a relatively rare situation) so the genotype caller
   treats it with appropriate caution.

Genotype calling, the stage that turns per-molecule counts into two alleles with a
confidence, was quite rough in the original ScaleHD (in hindsight). As we were discovering aspects about the
biology of HTT during development of the original ScaleHD, the goalposts for genotyping shifted often.
Also I was not as skilled as I am today. Now each candidate genotype is a mixture of two alleles, 
each smeared by a PCR slippage/stutter kernel whose shape depends on
repeat length (calibrated on the ScaleHD 1.x training matrix, so the model again may change after access to real data
or up-to-date scientific insight). Candidates are compared by how well they statistically explain every molecule, 
which gives a posterior probability for the call plus flags for the cases that deserve a manual look.

Read about the statistical model in [MODEL.md](MODEL.md).
