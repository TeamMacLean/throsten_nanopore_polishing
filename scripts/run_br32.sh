#!/bin/bash

cd /tsl/scratch/shrestha/thorsten_nanopre_assembly_polishing

bash scripts/pipeline.sh  BR32 /tsl/scratch/langnert/Minichromosomes/GEMO_project/Nanopore_reads/BR32_FG1801_02_unfiltered.fastq.gz /tsl/scratch/langnert/Minichromosomes/Minichromosomes/Genomes/Magnaporthe_nanopore_genomes_2018_Dec/BR32_genome_21.fasta /tsl/scratch/langnert/Minichromosomes/GEMO_project/Illumina_reads/BR32_ERR747964.fastq.gz
