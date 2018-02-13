.. _sect_outputhierarchy:

Output Hierarchy
================================

ScaleHD's output is fairly straightforward. For each job that you run, there will be one master output folder containing all data produced by the pipeline. If the user specified a '-j' argument in the command line interface, the output folder will be named as such; if not, the folder will be named ScaleHD-DateTime.

Within the master folder, there will be a few instance-wide report files:

 * AlignedDistributions.csv
 * InstanceGraphs.pdf
 * InstanceReport.csv
 * UtilisedConfiguration.xml

AlignedDistributions contains a matrix of read count distributions for each allele, in each sample processed by the pipeline. Each individual allele's distribution will be aligned to the same index within this matrix; the peak, or n-value, of each allele is at the same index. This allows the user an easy method by which to cross-compare levels of somatic mosaicism.

InstanceGraphs is a series of graphs rendered for each sample. These graphs plot the CCG distribution, and respective CAG distributions for each allele in a sample. This allows the user to easily search for samples that may require manual inspection and view read count distributions quickly.

InstanceReport is the 'main' report output from ScaleHD, and displays allele genotypes and relevant statistics, data characteristics and errors (where applicable). A detailed explanation of what each of the attributes/flags in this file explains, please see :ref:`sect_definitions`.

UtilisedConfiguration is a literal copy of the ScaleHD configuration file specified as input by the user when the job was launched.