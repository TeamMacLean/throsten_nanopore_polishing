## Feb 1

Obtained datasets for all samples. Studied the help page for running minimap2, racon and pilon.

I started to code a script for nanopore assembly polishing, one for each sample and a master script to run the pipeline.

Starting running 5 samples nanopore assembly polishing. 


## Feb 5 2019

Racon accepts PAF format file from minimap as shown in the help. Minimap2 outputs PAF format by default however, racon was not accepting the PAF file from minimap.

There is no way to change PAF to SAM format as PAF has missing CIGAR data.

I decided to rerun the job and this time minimap2 to produce sam format files to input in racon.

I made changes in the script and then sbatched the jobs to HPC. All jobs in pending at the moment.

## Feb 11 2019

Racon consensus generation for third round get error. Therefore, racon step is done for 2 rounds only. 

## Feb 13 2019

Racon and pilon assembly polishing completed for all samples.

## April 23 2019

Running racon and pilon on second version of assembly
