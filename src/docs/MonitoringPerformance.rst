.. _sect_perfmon:

Monitoring Performance
================================

If you're running ScaleHD in verbose mode (-v/--verbose), the pipeline will print to the terminal any updates on stages for each sample that is being processed. This involves a colour code, for success (green), warning/information (yellow) and failure (red). However, if you're more interested in progress for each individual sample as it is being processed, you will need to look at additional software. Personally, I find it easiest to use the program *htop*, which you can read about here: https://hisham.hm/htop/. Running htop will allow you to see what third party binaries that ScaleHD has launched (i.e. what stage of the sample processing it is on), as well as resource usage.

Normally, we just leave a job running overnight and return to it being finished (or almost finished) the next morning.