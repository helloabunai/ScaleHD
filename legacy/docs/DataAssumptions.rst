.. _sect_dataassume:

Data Assumptions
================================

ScaleHD makes a few assumptions about the input data which you provide it. Details about those assumptions can be found in this section, to aid the user in preparing their input into an appropriate format for ScaleHD to process it successfully.

File Input
~~~~~~~~~~~

ScaleHD utilises both forward and reverse reads from amplicon sequencing in order to genotype a sample effectively. In order to determine whether an individual file is a sample's forward or reverse data, ScaleHD looks for the substring "_R1" (forward) and "_R2" (reverse) at the end of each file name. Thus, user input should follow this behaviour within the specified input folder.

.. image:: img/file-name.png

In addition to the file name, ScaleHD expects any input folder to contain an even number of files, and only sequence input data (either *.fastq.gz, or unzipped) within said input folder. To be verbose, valid input file extensions are as follows:

 * sample_name.fastq.gz
 * sample_name.fq.gz
 * sample_name.fastq
 * sample_name.fq

No other files will be considered valid. If there is a non-even number of input files present (i.e. not every sample has two files, R1 and R2), ScaleHD will not run.

Multiplexed Data
~~~~~~~~~~~~~~~~

If your sequence data originates from a samplelibrary combined with other non-*HTT* data, then your data will need to be de-multiplexed before it can be used within ScaleHD. In order to demultiplex data, we typically use Cutadapt (http://cutadapt.readthedocs.io/en/stable/) to select *HTT* reads only on the basis of the presence of the PCR primer binding site in the reads. Even for non-multiplexed data (*HTT* only data) it is useful to use Cutadapt to collect reads that start with the PCR primer binding site. There are two benefits in doing this:

 * Removing the reads that, for some reason, would not start with a PCR primer binding site.
 * Trimming all the reads at the same position within the PCR primer binding site. This allows processing of reads which start at the same position, and also allows for trimming of the hetereogeneity spacer we often use for sequencing library preparation.

Due to colleague unfamiliarity with Bash scripting, it was easier to write a simple python wrapper to run Cutadapt in batch mode, for particular barcodes which our lab uses. See BatchAdapt at https://github.com/helloabunai/Batchadapt. Documentation for this small wrapper is in progress at the time of writing, so a simple instruction set is as follows:

 * Clone the Github repository
 * Run the setup script in a terminal: python setup.py install
 * Run "Batchadapt --help" to see all input options

More detailed explanation of the demultiplexing process is a work in progress and will be updated here when complete.

Reference Libraries
~~~~~~~~~~~~~~~~~~~

Due to the multiple reference design approach to determining legitimate alleles that ScaleHD uses, the pipeline expects an input reference library to conform to a certain standard.

.. _forward-reference:
For the *forward* reference library, there must be a total dimensional value of 4000 references within the library; CAG1-200 with CCG1-20. Each individual reference sequence label must also conform to a naming standard. For example, a reference sequence with CAG 17 and CCG 6 (i.e. CAG17_CCACAG1_CCGCCA1_CCG6_CCT2) must be labelled as:

::

  >17_1_1_6_2
  GCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACAGCCGCCACCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCCGCCGCAGGCACAGCCGCTGCT

.. _reverse-reference:
For the *reverse* reference library, we utilise a static 100CAG with varying CCG1-20. This is as we only utilise the reverse reads to determine CCG values and resultant zygosity. An example of an individual reverse reference is:

::

  >100_1_1_15_2
  CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACAGCCGCCACCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCCGCCGCAGGCACAGCCGCTGCT

Basically, the only thing that is able to be altered about a user's input reference library is the forward and reverse flanks of the sequence. This is an intentional design decision, and ScaleHD will not function properly with alternatively styled reference libraries. Please see my other software, RefGeneratr, at https://github.com/helloabunai/RefGeneratr. This package allows for generating custom reference libraries with ease, and is actually installed as a dependency for ScaleHD.

An example of valid reference libraries can be seen on ScaleHD's github page, under /src/ScaleHD/config/4k-HD-INTER.fas and /src/ScaleHD/config/20TypicalReverse.fasta. Valid file extensions are:

 * reference_name.fasta
 * reference_name.fas
 * reference_name.fa
