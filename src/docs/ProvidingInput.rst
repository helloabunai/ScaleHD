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

ScaleHD's main arguments, rather than filling out cumbersome details on a command line (where editing can be frustrating), instead uses an XML document. This allow allows for an easy mechanism by which to reproduce results of a processing job -- if you use the same configuration file, you can guarantee that data produced at the end will be identical.

Here is an example configuration file:

::

  <config data_dir="/path/to/fastq/" forward_reference="/path/to/fa" reverse_reference="/path/to/fa">
  <instance_flags quality_control="True" sequence_alignment="True" atypical_realignment="True" genotype_prediction="True" snp_calling="False"/>
	  <trim_flags trim_type="Adapter" quality_threshold="5" adapter_flag="-a" forward_adapter="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" reverse_adapter="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" error_tolerance="0.39"/>
  <alignment_flags min_seed_length="19" band_width="100" seed_length_extension="1.5" skip_seed_with_occurrence="500" chain_drop="0.50" seeded_chain_drop="0" seq_match_score="1" mismatch_penalty="4" indel_penalty="6,6" gap_extend_penalty="4,4" prime_clipping_penalty="5,5" unpaired_pairing_penalty="17"/>
  <prediction_flags plot_graphs="True"/>
  </config>

Within the **<config>** branch, there are three attributes to which the user must assign a value. *data_dir* must point to your input folder, consisting of an even number of input data files (see: :ref:`sect_dataassume`). *forward_reference* points to a *.fasta reference file, for which alignment is carried out on forward reads (see: :ref:forward_reference). *reverse_reference* points to a *.fasta reference file, for reverse alignment (see: :ref:reverse_reference).

**<instance_flags>** determines which stage(s) of ScaleHD that the user wishes to run. These are all simple True/False boolean options, where 'False' means a stage will not be processed. Currently, snp_calling is still in development and testing, and as such does not execute regardless of which boolean you input here. Other stages all function as expected.

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
+----------------------------------------------+
| prime_clipping_penalty    | -L [INT, INT]    |
+----------------------------------------------+
| unpaired_pairing_penalty  | -U <INT>         |
+----------------------------------------------+

**<prediction_flags>** is a remnant of previous development implementations, where we utilised different statistical approaches to genotyping and the end-user was expected to modify the parameters of a predictive model for improved accuracy. However, we have radically changed our design approach in this implementation of the pipeline and there is only one flag remaining, which decides whether you want graphs to be rendered or not. However, this has not been re-implemented within this approach to ScaleHD and as a result, doesn't do anything if you choose 'False'. ¯\_(ツ)_/¯.




