.. _sect_seqqc:

Sequence Quality Control (SeqQC)
================================

Assuming that all dependencies are functioning and input data is valid, the first stage of ScaleHD is Sequence Quality Control (SeqQC). This stage is utilised to attempt to provide further stages with the best possible quality data for aligning and genotyping.

As of version 0.312, we have included a pre-processing stage for demultiplexing sequence reads based on sequencing adapters. This is done via a wrapper for Cutadapt, called batchadapt -- it is a lazy way for my colleagues to easily run cutadapt for sequence trimming on a folder of input files without having to think about bash scripts. As such, batchadapt has been added to the requirements for this package, and is automatically installed from pip along with other dependencies.

Initially, the pipeline will execute FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on input data. This provides a visual report on raw sequence quality, which may inform users of basic characteristics of their sequencing data. 

The main focus of this stage utilises Cutadapt (http://cutadapt.readthedocs.io/en/stable/). In a ScaleHD configuration file, the user can specify three types of trimming to be carried out: quality, adapter or both. Quality trimming, if chosen, is recommended to be carried out before adapter trimming -- if the user has chosen to use both types of trimming via ScaleHD, this is the order in which commands will be executed. Quality trimming will remove low-quality reads from an input data, providing any further analysis stages with the highest possible quality data available. See http://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming for more information.

Adapter trimming takes a literal sequence specified by the user, and trims all target reads of this sequence. This stage can have quite a profound impact on the downstream quality of genotyping, as sequence alignment is affected by proper adapter trimming by quite a noteworthy margin. If you are going to trim adapters, please make sure you are using the correct sequence for your data, and the correct positioning!

.. note::
    Cutadapt adapter trimming requires additional reading to fully understand what action(s) is(are) executed. The application uses additional arguments to indicate what position on the trimming targets that any given adapter should be trimmed from (5 prime, 3 prime, anchored). Please see Cutadapt documentation on removing adapters for more information: http://cutadapt.readthedocs.io/en/stable/guide.html#removing-adapters.

Once complete, the trimmed data is then stored temporarily within the current sample's output folder. The input data is not moved or altered in any way; all further processing is carried out on the trimmed data, which is deleted upon sample completion. This avoids multiple copies of data being produced by each stage of the pipeline, and ensures that input data is not modified.