#!/bin/bash

source minimap2-2.13
source racon-1.2.0
cd /tsl/scratch/shrestha/thorsten_nanopre_assembly_polishing

sample=$1
nanopore_reads=$2
nanopore_assembly=$3
illumina_reads_1=$4
illumina_reads_2=$5


mkdir -p results/${sample}

if test ! -e results/${sample}/minimap2_nanopore_overlap1.sam.gz then
        echo "Running minimap2 round 1"
        cmd="minimap2 -a -x ava-ont -c -t 32 ${nanopore_assembly} ${nanopore_reads} | gzip - > results/${sample}/minimap2_nanopore_overlap1.sam.gz";  echo $cmd; $cmd
else
        echo "*** Minimap2 overlap round 1 completed ***"
fi

if test ! -e results/${sample}/racon_round1_output &&  test results/${sample}/minimap2_nanopore_overlap1.sam.gz -nt results/${sample}/racon_round1_output ; then
        echo "Running racon round 1"
        cmd="racon --threads 32 --include-unpolished  ${nanopore_reads} results/${sample}/minimap2_nanopore_overlap1.sam.gz ${nanopore_assembly} > results/${sample}/racon_round1_output";  echo $cmd; $cmd
else
        echo "*** Racon round 1 polishing completed ***"
fi

echo "*** generating overlap information after round 1 racon ***"

if test ! -e results/${sample}/minimap2_nanopore_overlap2.sam.gz && results/${sample}/racon_round1_output -nt results/${sample}/minimap2_nanopore_overlap2.sam.gz ; then
        echo Running minimap2 round 2
        cmd="minimap2 -a -x ava-ont -c -t 32 results/${sample}/racon_round1_output ${nanopore_reads} | gzip -  > results/${sample}/minimap2_nanopore_overlap2.sam.gz";  echo $cmd; $cmd
else
        echo "*** Minimap2 overlap round 2 completed ***"
fi

if test ! -e results/${sample}/racon_round2_output &&  test results/${sample}/racon_round1_output -nt results/${sample}/racon_round2_output; then
        echo Running racon round 2
        cmd="racon --threads 32 --include-unpolished ${nanopore_reads} results/${sample}/minimap2_nanopore_overlap2.sam.gz results/${sample}/racon_round1_output > results/${sample}/racon_round2_output";  echo $cmd; $cmd
else
        echo "*** Racon round 2  polishing completed ***"
fi


echo " *** generating overlap information after round 2 racon ***"

if test ! -e results/${sample}/minimap3_nanopore_overlap2.sam.gz && results/${sample}/racon_round2_output -nt results/${sample}/minimap3_nanopore_overlap2.sam.gz; then
        echo Running minimap2 round 3
        cmd="minimap2 -a -x ava-ont -c -t 32 results/${sample}/racon_round2_output ${nanopore_reads} | gzip - > results/${sample}/minimap2_nanopore_overlap3.sam.gz";  echo $cmd; $cmd
else
        echo "*** Minimap2 overlap round 3 completed ***"
fi

if test ! -e results/${sample}/racon_round3_output &&  test results/${sample}/racon_round2_output -nt results/${sample}/racon_round3_output; then
        echo Running racon round 3
        cmd="racon --threads 32 --include-unpolished ${nanopore_reads} results/${sample}/minimap2_nanopore_overlap3.sam.gz results/${sample}/racon_round2_output > results/${sample}/racon_round3_output";  echo $cmd; $cmd
else
        echo " **** Racon round 3 polishing completed ***"
fi


## next step map illumina reads to Nanopore polished assembly from above steps

source bowtie2-2.2.9
source pilon-1.22
source samtools-1.9

echo "*** running bowtie and pilon ***"

if test ! -e results/${sample}/racon_round3_output; then
        cmd=$(bowtie2-build results/${sample}/racon_round3_output results/${sample}/racon_round3_output);  echo $cmd; $cmd
fi

if test ! -e result/${sample}/illumina_aligned_to_assembly_sorted.bam; then
        if [ "$illumina_reads_2" == "" ]; then
                cmd="bowtie2 --threads 8 --qc-filter --no-unal -x results/${sample}/racon_round3_output -U $illumina_reads_1 | samtools view -b  | samtools sort -o results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam";  echo $cmd; $cmd
                cmd="pilon --genome results/${sample}/racon_round3_output --unpaired results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam --outdir results/${sample} --output ${sample}_pilcon_polished --changes --threads 16";  echo $cmd; $cmd
        else
                cmd="bowtie2 --threads 8 --qc-filter --no-unal -x results/${sample}/racon_round3_output -1 $illumina_reads_1 -2 $illumina_reads_2 | samtools view -b  | samtools sort -o results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam";  echo $cmd; $cmd
                cmd="pilon --genome results/${sample}/racon_round3_output --frags results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam --outdir results/${sample} --output ${sample}_pilcon_polished --changes --threads 16"; echo $cmd; $cmd
        fi
fi
