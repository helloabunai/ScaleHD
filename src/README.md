ScaleHD: Automated Huntington Disease genotyping
=========================================================
ScaleHD is a package for automating the process of genotyping microsatellite repeats in Huntington Disease data.
We utilise machine learning approaches to take into account natural data 'artefacts', such as PCR slippage and somatic
mosaicism, when processing data. This provides the end-user with a simple to use platform which can robustly predict genotypes from input data.

By default, input is a pair of unaligned .fastq sequence data -- both forward and reverse reads, per sample. We utilise both forward and reverse
reads in order to reduce the complex dimensionality issue posed by Huntington Disease's multiple repeat tract genetic structure. Reverse reads allow
us to determine the current sample's CCG state -- this provides us with a mechanism by which to more easily call the entire genotype. Forward reads
are utilised in a similar approach, to determine the CAG and intervening structure.

The general overview of the application is as follows:
1) Input FastQ files are subsampled, if an overwhelming number of reads are present. This can be overruled with the -b flag.
2) Sequence quality control is carried out per the user's instructions. We reccomend triming of any 5-prime spacer+primer combinations, for optimal alignment.
3) Alignment of these files, to a typical HD structure (CAG_1_1_CCG_2) reference, is carried out.
4) Assemblies are scanned with Digital Signal Processing to detect any possible atypical structures (e.g. CAG_2_1_CCG_3).
4.1) If no atypical alleles are detected, proceed as normal.
4.2) If atypical alleles are detected, a custom reference is generated, and re-alignment to this is carried out.
5) With the appropriate allele information and sequence assembly(ies) present, sampled are genotyped.
6) Output is written for the current sample; the procedure is repeated for the next sample in the queue (if present).


What's New
==========
* Added an n-aligned matrix of repeat-count distributions, on a SHD instance-wide basis.
* Instances of SHD will output the utilised configuration file with other results.
* Removed the -b/--boost flag, and made subsampling the default behaviour (given acceptable read-count in raw input data)
* Added the -b/--broadscope flag, which forces alignment and DSP to be executed on all present reads (i.e. no subsampling).
* Added the -e/--enshrine flag, which forces SHD to retain all aligned reads which are not uniquely mapped (which are removed by default).
* Implemented DSP to function within the intervening sequence, rather than utilising a string derision method.
* Added many report flags for SHD's instance report output -- more on this in the Output section of this readme.


Installation Prerequisites
==========================
If you do not have sudo access (to install requisite packages), you should run ScaleHD in a user-bound local python environment,
 or discrete installation. This guide will assume you have sudo access. However, we detail an extra stage on setting up a local
 python environment for use with ScaleHD.

0. (Optional 1 - no sudo) Python 2.7 Setup
    ~~~~
    $ cd desired-directory
    $ tar jvzf Python-2.7.tar.bz2
    $ cd Python-2.7
    $ ./configure --enable-shared --prefix=/your/custom/installation/path
    $ make
    $ make install
    ~~~~

0. (Optional 2 - no sudo) Bash profile edit.. in your ~/.bash_profile file
    ~~~~
    $ export PATH=/your/custom/installation/path/bin:$PATH
    $ export LD_LIBRARY_PATH=/your/custom/installation/path/lib:$LD_LIBRARY_PATH
    ~~~~

1. Get PIP if not already installed!
    ~~~~
    $ wget https://bootstrap.pypa.io/get-pip.py
    $ python ~/path/to/get-pip.py
    ~~~~

2. Install Cython/Scipy stack separately (Setuptools seems to install incorrectly..)
    ~~~~
    $ pip install cython
    $ pip install scipy
    $ pip install numpy
    ~~~~

3. Install ScaleHD from src (pip coming soon...)
    ~~~~
    $ cd ~/path/to/ScaleHD/src/
    $ python setup.py install
    ~~~~

4. Install required third-party binaries. Please make sure any binaries you do install are included on your $PATH so that they can be found by your system. 
**Please note**, ScaleHD will utilise GNU WHICH/TYPE to determine if a command is on your $PATH. If either WHICH/TYPE or a dependency is missing, ScaleHD will inform you and exit.
    ~~~~
    Quality Control:
        Cutadapt
        FastQC (Java required)
    Alignment:
        BWA        
        SeqTK
        Samtools
        Generatr (setup.py will install this for you)
    Genotyping:
        R        
        Samtools
        Generatr (as above)
        Picard (alias required*)
        GATK (alias required*)
    ~~~~
*aliases are required for certain third party software which are not distributed as installable binaries. An example of an alias would
look like:

    alias gatk="java -jar /Users/home_dir/Documents/Builds/GenomeAnalysisTK.jar"

5. Check that libxml2-dev and libxslt-dev are installed...

Usage
=====
General usage is as follows:

    $ scalehd [-h/--help] [-v] [-c CONFIG] [-t THREADS] [-e] [-b] [-g] [-j "jobname"] [-o OUTPUT]
    e.g.
    $ scalehd -v -c ~/path/to/config.xml -t 12 -j "ExampleJobPrefix" -o ~/path/to/master/output

ScaleHD flags are:

    -h/--help:: Simple help message explaining flags in detail
    -v/--verbose:: Enables verbose mode in the terminal (i.e. shows user feedback)
    -c/--config:: Will execute all settings specified in the given ArgumentConfig.xml [filepath].
    -t/--threads:: Number of threads to utilise. Mainly will affect alignment performance [integer].
    -e/--enshrine:: Forces aligned reads which are not uniquely mapped to be retained; default behaviour without this flag removes said reads.
    -b/--broadscope:: Forces subsampling of raw and aligned reads to be disabled.
    -g/--groupsam:: Groups all aligned assemblies generated into one output folder, with appropriate sample names. If not specified, assemblies will be left in the sample's specific output subfolder.
    -j/--jobname:: Specifies a prefix to use for the root output directory. Optional. If you specify a JobName that already
    exists within your specified -o output folder, ScaleHD will prompt the user to decide if they wish to delete the pre-existing folder and replace.
    -o/--output:: Desired output directory.

Data Primer
===========
A short note on the requirements of filenames/structure for ScaleHD to function. A sample's filename (here, ExampleSampleName) must adhere to the following structure:
    
    ExampleSampleName_R1.fastq
    ExampleSampleName_R2.fastq
    
You must utilise both forward (R1) and reverse (R2) reads, per sample pair. If the respective files do not end in _R1.fastq (.fq) or _R2.fastq (.fq), ScaleHD will not run correctly.
Since this is a highly HD specific application, we can offer some insight into providing the best approaches for genotyping. Due to the similarity of both repeat tracts in HD (CAG and CCG), which
are flanking an intervening sequence, that in itself is highly similar to both regions, alignment can be fussy about your input data. Thus, we highly recommend trimming any spacers or primers present
on the 5Prime end of your reads; this enables reads to start at the same position and provides the aligner with a more discrete boundary between the different HD repeat tracts.

Individual settings for different stages in ScaleHD are set within a configuration XML document. The particular acceptable data types/ranges for each parameter varies. The configuration XML document for ScaleHD settings must also adhere to the following structure:

    <config data_dir="/path/to/reads/" forward_reference="/path/to/forward/ref_seq.fa" reverse_reference="/path/to/reverse/ref_seq.fa">
      <instance_flags quality_control="BOOL" sequence_alignment="BOOL" atypical_realignment="BOOL" genotype_prediction="BOOL", snp_calling="BOOL"/>
      <trim_flags trim_type="x" quality_threshold="x" adapter_flag="x" adapter="x" error_tolerance="x"/>
      <alignment_flags min_seed_length="x" band_width="x" seed_length_extension="x" skip_seed_with_occurrence="x" chain_drop="x" 
      seeded_chain_drop="x" seq_match_score="x" mismatch_penalty="x" indel_penalty="x" gap_extend_penalty="x" prime_clipping_penalty="x" 
      unpaired_pairing_penalty="x"/>
      <prediction_flags plot_graphs="BOOL"/>
    </config>
    
With each parameter data type/rule being as follows:
    
    CONFIG
        data_dir: Must be a real path, with an even number of ONLY *.fastq or *.fq files within.
        forward_reference: Must be a real reference file (*.fasta, *.fa or *.fas).
        reverse_reference: See forward_reference.
    INSTANCE
        quality_control: Boolean, TRUE/FALSE
        sequence_alignment: Boolean, TRUE/FALSE
        atypical_realignment: Boolean, TRUE/FALSE
        genotype_prediction: Boolean, TRUE/FALSE
        snp_calling: Boolean, TRUE/FALSE
    TRIM
        trim_type: String, "Quality", "Adapter" or "Both"
        quality_threshold: Integer, within the range 0-38
        adapter_flag: String, one of: '-a','-g','-a$','-g^','-b'. ([See Cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html#removing-adapters))
        adapter: String, consisting of only 'A','T','G','C'
        error_tolerance: Float, within the range of 0.0 to 1.0 (in 0.01 increments).
    ALIGNMENT
        All flags present are direct equivalents of parameters present in BWA-MEM. 
        See [the BWA manual for more information](http://bio-bwa.sourceforge.net/bwa.shtml).
    PREDICTION
        plot_graphs: Boolean, TRUE/FALSE
    
Output
======
A brief overview of flags provided in the output is as follows:

    SampleName:: The extracted filename of the sample that was processed.
    Primary/Secondary GTYPE:: Allele genotype in the format CAG_x_y_CCG_z
    Status:: Atypical or Typical structure
    BSlippage:: Slippage ratio of allele's read peak ('N minus 2' to 'N minus 1)', over 'N'.
    Somatic Mosaicism:: Mosaicism ratio of allele's read peak ('N plus 1' to 'N plus 10'), over 'N'
    Confidence:: Confidence in genotype prediction (0-100).
    Exception Raised:: If, during a particular stage of the pipeline, exceptions caused the processing to fail, this flag will inform the user in which stage it crashed.
    Homozygous Haplotype:: If True, both alleles have an identical genotype.
    Neighbouring Peaks:: If True, both alleles exist within the same CCG distribution, neighbouring each other.
    Diminished Peaks:: If True, an expanded peak has very few reads and was detected independently. Manual inspection recommended.
    Novel Atypical:: If True, an intervening sequence structure that has not been readily observed before was detected. Manual inspection recommended.
    Alignment Warning:: If True, determining the CCG value(s) returned more peaks than is 'possible'. Manual inspection recommended.
    Atypical Alignment Warning:: In the case of atypical re-alignment, particularly awful alignment quality can return more than one peak; which should not happen.
    CCG Rewritten:: CCG was rewritten from the FOD-derived value -- i.e. DSP overwrote the FOD results.
    CCG Zygosity Rewritten:: A sample (aligned to a typical reference) that was heterozygous (CCG), was detected to be an atypical homozygous (CCG) sample.
    CCT Uncertainty:: The most common CCT 'sizes' returned from DSP were too similar in count (e.g. CCT2 == 54%, CCT3 == 46%) to be certain.
    SVM Failure:: SVM CCG zygosity calling was incorrect, as a result of the resultant confusion matrix providing differing results from a brute force ratio check. Manual inspection highly recommended.
    Differential Confusion:: The allele sorting algorithm is confused between a potential neighbouring peak, and a homozygous haplotype. Manual inspection highly recommended.
    Peak Inspection Warning:: At least one allele failed inspection on the repeat-count distribution the genotype(s) was(were) derived from. Common in very low read count samples/poor sequencing.
    Low Distribution Reads:: A warning which is triggered when at least one allele's repeat count distribution contains an unappealingly low number of reads.
    Low Peak Reads:: A fatal warning which is triggered when, in a given repeat count distribution, the returned N value contains a very low number of reads. Manual inspection highly recommended.