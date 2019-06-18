#!/bin/bash

#bash scripts/pipeline.sh FR13 /tsl/scratch/langnert/Minichromosomes/GEMO_project/Nanopore_reads/FR13_FG1801_05_unfiltered.fastq.gz /tsl/scratch/shrestha/thorsten_nanopre_assembly_polishing/data/FR13/FR13_32.fa /tsl/scratch/langnert/Minichromosomes/GEMO_project/Illumina_reads/FR13_APS_EOSU_2_62FJ7AAXX.fastq.gz

source ruby-2.3.1
# read data 
R1=/tsl/scratch/langnert/Minichromosomes/GEMO_project/Illumina_reads/FR13_APS_EOSU_2_62FJ7AAXX.fastq.gz
nanopore=/tsl/scratch/langnert/Minichromosomes/GEMO_project/Nanopore_reads/FR13_FG1801_05_unfiltered.fastq.gz
assembly=/tsl/scratch/shrestha/thorsten_nanopre_assembly_polishing/data/FR13/FR13_32.fa

rake -f scripts/pipeline.rake sample=FR13 nanopore_reads=${nanopore} assembly=${assembly} R1=$R1
