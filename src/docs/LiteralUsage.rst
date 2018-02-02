.. _sect_literalusage:

Running ScaleHD
================================

An example of running the program in a "default" state (max threads, remove non-unique reads, subsample high read count) would look like the following:

::
  ScaleHD -v -c ~/Documents/Work/ScaleHDTests/ArgumentConfig.xml -j "hello_documentation" -o ~/Documents/Work/ScaleHDTests/Output

Once this command is entered into a terminal, ScaleHD will read the input from ArgumentConfig.xml, process the input data specified there, and genotype all samples where possible. This is all you need to do to run ScaleHD; when the input is set up correctly, no input is required from the end user.