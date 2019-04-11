.. _sect_input:

Specifying User Input(s)
================================

ScaleHD is a Command Line Interface (CLI) based application, and users will interact with it via a terminal environment. Once installed, the program can be launched by simply issuing the command 'scalehd' in a terminal. This, by default, will not do anything, but will list all the available input options to the user. Here we describe what these options are, how the user can modify them and what those changes will do to the functionality of ScaleHD.

+----------------+---------------+---------------------------------------------------------+
| Short argument | Long argument | Function                                                |
+================+===============+=========================================================+
| -h             | --help        | Prints a help message                                   |
+----------------+---------------+---------------------------------------------------------+
| -v             | --verbose     | Enables terminal feedback for the user (def: blank)     |
+----------------+---------------+---------------------------------------------------------+
| -c {string}    | --config      | Specify a directory to your configuration file.         |
+----------------+---------------+---------------------------------------------------------+
| -t {integer}   | --threads     | Number of CPU threads to use (def: system max).         |
+----------------+---------------+---------------------------------------------------------+
| -p             | --purge       | Remove all output files except the HTML5-based report.  |
+----------------+---------------+---------------------------------------------------------+
| -s             | --simple      | Enable simple 95% confidence interval genotyping output.|
+----------------+---------------+---------------------------------------------------------+
| -e             | --enshrine    | Do not remove non-uniquely mapped read from SAM files.  |
+----------------+---------------+---------------------------------------------------------+
| -b             | --broadscope  | Do not subsample fastq data with high read counts.      |
+----------------+---------------+---------------------------------------------------------+
| -g             | --groupsam    | Output all sorted SAM files into one job-wide directory.|
+----------------+---------------+---------------------------------------------------------+
| -j {string}    | --jobname     | Customises the output folder name for this job.         |
+----------------+---------------+---------------------------------------------------------+
| -o {string}    | --output      | Specify a directory you wish output to be directed at   |
+----------------+---------------+---------------------------------------------------------+

Configuration File
~~~~~~~~~~~~~~~~~~

ScaleHD's main arguments, rather than filling out cumbersome details on a command line (where editing can be frustrating), instead uses an XML document. This allows for an easy mechanism by which to reproduce results of a processing job -- if you use the same configuration file, you can guarantee that data produced at the end will be identical.

Here is an example configuration file:

::

<config data_dir="~/path/" forward_reference="~/path/" reverse_reference="~/path/">
	<instance_flags demultiplex="False" quality_control="True" sequence_alignment="True" atypical_realignment="True" genotype_prediction="True" snp_calling="True"/>
	<demultiplex_flags forward_adapter="GCGACCCTGG" forward_position="5P" reverse_adapter="GCAGCGGCTG" reverse_position="5P" error_rate="0" min_overlap="10" min_length="" max_length=""/>
	<trim_flags trim_type="Adapter" quality_threshold="5" adapter_flag="-a" forward_adapter="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" reverse_adapter="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" error_tolerance="0.39"/>
	<alignment_flags min_seed_length="19" band_width="100" seed_length_extension="1.5" skip_seed_with_occurrence="500" chain_drop="0.50" seeded_chain_drop="0" seq_match_score="1" mismatch_penalty="4" indel_penalty="6,6" gap_extend_penalty="6,6" prime_clipping_penalty="5,5" unpaired_pairing_penalty="17"/>
	<prediction_flags snp_observation_threshold="2" quality_cutoff="0"/>
</config>


Within the **<config>** branch, there are three attributes to which the user must assign a value. *data_dir* must point to your input folder, consisting of an even number of input data files (see: :ref:`sect_dataassume`). *forward_reference* points to a .fasta reference file, for which alignment is carried out on forward reads. *reverse_reference* points to a .fasta reference file, for reverse alignment.

**<instance_flags>** determines which stage(s) of ScaleHD that the user wishes to run. These are all simple True/False boolean options, where 'False' means a stage will not be processed. Other stages all function as expected.

**<demultiplex_flags>** provides input arguments for demultiplexing, if chosen. These are arguments for Batchadapt, which is a simple wrapper I wrote for a colleague to run Cutadapt (http://cutadapt.readthedocs.io/en/stable/) on a folder of input sample pairs. x_position flags must be '5P', '3P' or 'AP'. Adapter strings must only contain A,T,G,C. Other arguments are positive integers, only.

**<trim_flags>** is for specifying options that get passed to cutadapt for sequence adapter/quality trimming. For an explanation of these flags, you can find the cutadapt documentation at: http://cutadapt.readthedocs.io/en/stable/.

.. note::
    *forward_adapter* and *reverse_adapter* are quite important, and you should carefully select how much adapter you wish to trim from your sequences. Incorrect trimming can have a notable downstream effect on the quality of resultant alignments, and thus, genotypes.

**<alignment_flags>** is for specifying options that get passed to bwa-mem for sequence alignment. A detailed explanation of all bwa-mem flags can be found in the page's documentation at: http://bio-bwa.sourceforge.net/bwa.shtml. However, you may be confused about which command line flags correspond to which XML fields in ScaleHD's input configuration. An explanation:

+---------------------------+------------------+
| ScaleHD XML               | BWA-MEM argument |
+===========================+==================+
| min_seed_length           | -k <INT>         |
+---------------------------+------------------+
| band_width                | -w <INT>         |
+---------------------------+------------------+
| seed_length_extension     | -r <FLOAT        |
+---------------------------+------------------+
| skip_seed_with_occurrence | -c <INT>         |
+---------------------------+------------------+
| chain_drop                | -D <FLOAT>       |
+---------------------------+------------------+
| seeded_chain_drop         | -W <INT>         |
+---------------------------+------------------+
| seq_match_score           | -A <INT>         |
+---------------------------+------------------+
| mismatch_penalty          | -B <INT>         |
+---------------------------+------------------+
| indel_penalty             | -O [INT, INT]    |
+---------------------------+------------------+
| gap_extend_penalty        | -E [INT, INT]    |
+---------------------------+------------------+
| prime_clipping_penalty    | -L [INT, INT]    |
+---------------------------+------------------+
| unpaired_pairing_penalty  | -U <INT>         |
+---------------------------+------------------+

**<prediction_flags>** now has a function! wow! Since SNP calling has been implemented, the user has the option of choosing the cutoff for filtering what we consider to be a valid SNP. Values accepted are range(1,5). 1 being the most lenient value in determining a SNP as valid, 5 being the most harsh. I typically use 2. Do what you want.

Quality cutoff indicates a positive integer which is utilised as a further filter for variant validation.

The VCF "QUAL" is described in the latest specification as:

::

	QUAL - quality: Phred-scaled quality score for the assertion made in ALT. i.e. −10log10 prob(call in ALT is
    wrong). If ALT is ‘.’ (no variant) then this is −10log10 prob(variant), and if ALT is not ‘.’ this is −10log10
    prob(no variant). If unknown, the missing value should be specified. (Numeric)

The user can control the Phred-scaled score lower limit for a valid SNP to be included with ScaleHD output via this option. This is entirely up to the user, and how sensitive you want SNP calling to be.
If ScaleHD misses a mutation because specified quality cutoffs were too strict, that's not really ScaleHD's fault.
