#!/bin/bash

#bash scripts/pipeline.sh CD156 /tsl/scratch/langnert/Minichromosomes/GEMO_project/Nanopore_reads/CD156_FG1801_08_pass.fastq.gz /tsl/scratch/langnert/Minichromosomes/GEMO_project/Nanopore_polishing_version2/CD156/CD156_joined_28.fasta /tsl/scratch/langnert/Minichromosomes/GEMO_project/454_reads/ERR747968.fastq /tsl/scratch/langnert/Minichromosomes/GEMO_project/454_reads/ERR747969.fastq /tsl/scratch/langnert/Minichromosomes/GEMO_project/454_reads/ERR747970.fastq /tsl/scratch/langnert/Minichromosomes/GEMO_project/454_reads/ERR747971.fastq

source ruby-2.3.1

nanopore=/tsl/scratch/langnert/Minichromosomes/GEMO_project/Nanopore_reads/CD156_FG1801_08_pass.fastq.gz
assembly=/tsl/scratch/langnert/Minichromosomes/GEMO_project/Nanopore_polishing_version2/CD156/CD156_joined_28.fasta

# first read set
R1=/tsl/scratch/langnert/Minichromosomes/GEMO_project/454_reads/ERR747968.fastq
# second read set or read paired to R1
/tsl/scratch/langnert/Minichromosomes/GEMO_project/454_reads/ERR747969.fastq
# third read set, if you have or forward read of second paired end reads
R3=/tsl/scratch/langnert/Minichromosomes/GEMO_project/454_reads/ERR747970.fastq
# forth read set, if you have or reverse read to second paired end reads
R4=/tsl/scratch/langnert/Minichromosomes/GEMO_project/454_reads/ERR747971.fastq

rake -f scripts/pipeline.rake sample=CD156 nanopore_reads=/tsl/scratch/langnert/Minichromosomes/GEMO_project/Nanopore_reads/CD156_FG1801_08_pass.fastq.gz assembly=/tsl/scratch/langnert/Minichromosomes/GEMO_project/Nanopore_polishing_version2/CD156/CD156_joined_28.fasta R1=/tsl/scratch/langnert/Minichromosomes/GEMO_project/454_reads/ERR747968.fastq R2=/tsl/scratch/langnert/Minichromosomes/GEMO_project/454_reads/ERR747969.fastq R3=/tsl/scratch/langnert/Minichromosomes/GEMO_project/454_reads/ERR747970.fastq R4=/tsl/scratch/langnert/Minichromosomes/GEMO_project/454_reads/ERR747971.fastq CD156

rake -f scripts/pipeline.rake sample=CD156 nanopore_reads=${nanopore} assembly=${assembly} R1=${R1} R2=${R2} R3=${R3} R4=${R4} CD156

