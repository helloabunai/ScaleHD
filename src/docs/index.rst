ScaleHD - Automated Huntington Disease genotyping pipeline
========================

ScaleHD is a bioinformatics pipeline for use in broad-scope automated genotyping of Huntington Disease (HD) data. Existing software which aim to profile nucleotide repeat diseases (such as HD) typically function on a generalised level, not focusing on any one particular disease. In the case of HD, this can result in inaccuracies occurring in any genotypes produced, due to a plethora of biological phenomena which are present, complicating matters when it comes to genotyping HD. This software takes a direct approach to resolving these problems, and doing so in an unsupervised, automated manner.

ScaleHD takes a configuration XML document as input, which contains all required information for the instance of ScaleHD to run to completion. The details of this XML input are specified in :ref:`sect_input`. Raw data should be in FastQ format; both forward (_R1) and reverse (_R2) sequence reads are required for ScaleHD to function. ScaleHD will perform quality control, sequence alignment and genotyping on all read-pairs presented by the user as input. If a sample fails to produce a genotype at any given stage, for whatever reason, a debug log is created so the user can (hopefully) understand why.

The source for ScaleHD is open source, and is available on Github at http://www.github.com/helloabunai/ScaleHD/ 

Currently, we are on Version 0.252. While the base algorithm has been implemented, much improvement still remains. As such, there will be continual development of ScaleHD for the forseeable future. For more information on version changes, please check :ref:`developer-documentation`.

The documentation for this software is organised into the relevant sections.

 * :ref:`noteworthy-features`
 * :ref:`prerequisites-and-assumptions`
 * :ref:`using-pipeline`
 * :ref:`understanding-output`

If you want to contact me, I can be reached at alastair.maxwell@glasgow.ac.uk.

Development is undertaken at the University of Glasgow, Scotland, and is funded by the Cure Huntington Disease Foundation (CHDI): http://www.chdifoundation.org/

.. _noteworthy-features:

.. toctree::
   :maxdepth: 2
   :caption: Noteworthy Features

   Overview
   SeqQC
   Alignment
   AutomatedGenotyping

.. _prerequisites-and-assumptions:

.. toctree::
   :maxdepth: 2
   :caption: Prerequisites and Assumptions

   RequiredPackages
   HardwareRequirements
   DataAssumptions

.. _using-pipeline:

.. toctree::
   :maxdepth: 2
   :caption: Using the Pipeline

   ProvidingInput
   LiteralUsage
   MonitoringPerformance

.. _understanding-output:

.. toctree::
   :maxdepth: 2
   :caption: Understanding ScaleHD Output

   OutputHierarchy
   Definitions

.. _developer-documentation:

.. toctree::
   :maxdepth: 2
   :caption: Developer Documentation

   VersionChangelog
   PlannedFeatures

Thanks for using ScaleHD.