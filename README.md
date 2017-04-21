ScaleHD: Automated Huntington Disease genotyping
=========================================================
ScaleHD is a package for automating the process of genotyping microsatellite repeats in Huntington Disease data.
We utilise machine learning approaches to take into account natural data 'artefacts', such as PCR slippage and somatic
mosaicism, when processing data. This provides the end-user with a simple to use platform which can robustly predict genotypes from input data.

By default, input is an aligned .sam file (either through stdin, or user specified files/directories); only genotyping is carried out.
However, if you wish to use ScaleHD as a pipeline for unaligned reads, providing the software with a configuration XML file will allow for
quality control (trimming, demultiplexing) of raw reads, alignment, and then genotyping.

The general overview of the application (assuming use of all stages) is as follows:
1) Input FastQ files are subsampled, if specified. Reads are then treated for quality (trimming, scoring), given the user's parameters.
2) FastQC is carried out on the treated files, with reports available in a given sample's output folder.
3) Alignment of these files, to a typical HD structure (CAG_1_1_CCG_2) reference, is carried out.
4) Assemblies are scanned with Digital Signal Processing to detect any possible atypical structures (e.g. CAG_2_1_CCG_3).
4.1) If no atypical alleles are detected, proceed as normal.
4.2) If atypical alleles are detected, a custom tailored reference is generated, and re-alignment to this is carried out.
5) With the appropriate allele information and sequence assembly(ies) present, samples are genotyped.
6) Output is written for the current sample; the procedure is repeated for the next sample in the queue (if present).

What's New
==========
* Added -s/--subsample flag for subsampling of input FastQ sequences (random subsampling based on user-provided float, e.g. 0.2 = 20%)
* Added -b/--boost flag to increase DSP scanning by subsampling aligned assemblies (generic SAM subsampling given the following rules):
-- *If File_ReadCount > 20000: subsample = 25%.*
-- *Elif 20000 > File_ReadCount > 10000: subsample = 20%*
-- *Else subsample = 10%*
* Added -g/--groupsam flag, to output all generated alignment files into a dedicated assembly folder; as opposed to sample specific subfolders.
* Added -j/--jobname flag for customised prefix of root output directories.
* Added -p/--purgesam flag for removing all reads from an alignment file (generated or specified) which are not uniquely mapped.
* Restructured the entire alignment/genotyping process to utilise object-based allele structure to allow dynamic per-allele remapping for single atypical alleles.
* Transitioned from a Bowtie2 wrapper to a BWA-MEM wrapper in -c mode; XML/DTD tags updated to reflect this.

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
**Please note**, ScaleHD will utilise GNU TYPE to determine if a command is on your $PATH. If either TYPE or a dependency is missing, ScaleHD will inform you and exit.
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
    ~~~~

Usage
=====

General usage is as follows:

    $ scalehd [-h/--help] [-v] [-c CONFIG] [-t THREADS] [-p] [-g] [-s FLOAT] [-b] [-j "jobname"] [-o OUTPUT]
    e.g.
    $ scalehd -v -c ~/path/to/config.xml -t 12 -p -s 0.5 -j "ExampleJobPrefix" -o ~/path/to/master/output

ScaleHD flags are:

    -h/--help:: Simple help message explaining flags in detail
    -v/--verbose:: Enables verbose mode in the terminal (i.e. shows user feedback)
    -c/--config:: Will execute all settings specified in the given ArgumentConfig.xml [filepath].
    -t/--threads:: Number of threads to utilise. Mainly will affect alignment performance [integer].
    -p/--purgesam:: Enables the purging of reads which are not uniquely mapped to a reference. Optional.
    -g/--groupsam:: Groups all aligned assemblies generated into one output folder, with appropriate sample names. If not specified, assemblies will be left in the sample's specific output subfolder.
    -s/--subsample:: subsamples input FASTQ (unaligned) data by the specified float, as percentage [float].
    -b/--boost:: generic subsampling of SAM (aligned) data by a percentage, dependant on how many reads aligned.
    -j/--jobname:: Specifies a prefix to use for the root output directory. Optional. If you specify a JobName that already
    exists within your specified -o output folder, ScaleHD will prompt the user to decide if they wish to delete the pre-existing folder and replace.
    -o/--output:: Desired output directory.
