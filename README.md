ScaleHD: Automated Huntington Disease genotyping
=========================================================
ScaleHD is a package for automating the process of genotyping microsatellite repeats in Huntington Disease data.
We utilise machine learning approaches to take into account natural data 'artefacts', such as PCR slippage and somatic
mosaicism, when processing data. This provides the end-user with a simple to use platform which can robustly predict genotypes
from input data.

By default, input is an aligned .sam file (either through stdin, or user specified files/directories); only genotyping is carried out.
However, if you wish to use ScaleHD as a pipeline for unaligned reads, providing the software with a configuration XML file will allow for
quality control (trimming, demultiplexing) of raw reads, alignment, and then genotyping.

What's New
==========
Added -j/--jobname flag for customised prefix of root output directories.

Added -p/--purgesam flag for removing all reads from an alignment file (generated or specified) which are not uniquely mapped.

Transitioned from a Bowtie2 wrapper to a BWA-MEM wrapper in -c mode; XML/DTD tags updated to reflect this.

Installation Prerequisites
==========================

If you do not have sudo access (to install requisite packages), you should run ScaleHD in a user-bound local python environment,
 or discrete installation. This guide will assume you have sudo access. However, we detail an extra stage on setting up a local
 python environment for use with ScaleHD.

0. (Optional 1) Python 2.7 Setup
    ~~~~
    $ cd desired-directory
    $ tar jvzf Python-2.7.tar.bz2
    $ cd Python-2.7
    $ ./configure --enable-shared --prefix=/your/custom/installation/path
    $ make
    $ make install
    ~~~~

0. (Optional 2) Bash profile edit.. in your ~/.bash_profile file
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

4. Install required third-party binaries. Not all are required; if you only wish to genotype you only are required to have Samtools installed.
Please make sure any binaries you do install are included on your $PATH so that they can be found by your system.
    ~~~~
    Quality Control:
        FastQC (Java required)
        Cutadapt
    Alignment:
        BWA
        Samtools
    Genotyping:
        Samtools
    ~~~~

Usage
=====

General usage is as follows:

    ~~~~
    $ scalehd [-h/--help] [-v] (-b BATCH | -c CONFIG) [-t THREADS] [-p] [-j JOBNAME] [-o OUTPUT]
    i.e.
    $ scalehd -v -b ~/path/to/samfiles -t 12 -p -j "ExampleJobPrefix" -o ~/path/to/output
    ~~~~

ScaleHD flags are:

    ~~~~
    -h/--help:: Simple help message explaining flags in detail
    -v/--verbose:: Enables verbose mode in the terminal (i.e. shows user feedback)
    -b/--batch:: Batch mode. Genotyping on a batch of pre-aligned SAM files [directory].
    -c/--config:: Config mode. Will execute all settings specified in the given ArgumentConfig.xml [filepath].
    -t/--threads:: Number of threads to utilise. Mainly will affect alignment performance [integer].
    -p/--purgesam:: Enables the purging of reads which are not uniquely mapped to a reference. Optional.
    -j/--jobname:: Specifies a prefix to use for the root output directory. Optional.
    -o/--output:: Desired output directory.
    ~~~~
