#!/bin/bash
# usage : bash pipeline.sh samplename nanopore_reads assembly R1 R2

source minimap2-2.13
source racon-1.2.0
cd /tsl/scratch/shrestha/thorsten_nanopre_assembly_polishing

sample=$1
nanopore_reads=$2
nanopore_assembly=$3
illumina_reads_1=$4
illumina_reads_2=$5
R3=$6
R4=$7

mkdir -p results/${sample}

if test ! -e results/${sample}/minimap2_nanopore_overlap1.sam.gz; then
        echo "Running minimap2 round 1"
        minimap2 -a -x ava-ont -t 32 ${nanopore_assembly} ${nanopore_reads} | gzip > results/${sample}/minimap2_nanopore_overlap1.sam.gz
		echo "*** Minimap2 overlap round 1 completed ***"
else
        echo "*** Minimap2 overlap round 1 has already completed ***"
fi

if test ! -e results/${sample}/racon_round1_output.fasta &&  test results/${sample}/minimap2_nanopore_overlap1.sam.gz -nt results/${sample}/racon_round1_output.fasta ; then
        echo "Running racon round 1"
        racon --threads 32 --include-unpolished  ${nanopore_reads} results/${sample}/minimap2_nanopore_overlap1.sam.gz ${nanopore_assembly} > results/${sample}/racon_round1_output.fasta
		echo "*** Racon round 1 polishing completed ***"
else
        echo "*** Racon round 1 polishing has already completed ***"
fi

echo "*** generating overlap information after round 1 racon ***"

if test ! -e results/${sample}/minimap2_nanopore_overlap2.sam.gz && test results/${sample}/racon_round1_output.fasta -nt results/${sample}/minimap2_nanopore_overlap2.sam.gz; then
        echo Running minimap2 round 2
        minimap2 -a -x ava-ont -t 32 results/${sample}/racon_round1_output.fasta ${nanopore_reads}  | gzip > results/${sample}/minimap2_nanopore_overlap2.sam.gz
		echo "*** Minimap2 overlap round 2 completed ***"
else
        echo "*** Minimap2 overlap round 2 has already completed ***"
fi

if test ! -e results/${sample}/racon_round2_output.fasta &&  test results/${sample}/minimap2_nanopore_overlap2.sam.gz -nt results/${sample}/racon_round2_output.fasta; then
        echo Running racon round 2
        racon --threads 32 --include-unpolished ${nanopore_reads} results/${sample}/minimap2_nanopore_overlap2.sam.gz results/${sample}/racon_round1_output.fasta > results/${sample}/racon_round2_output.fasta
		echo "*** Racon round 2 polishing completed ***"
else
        echo "*** Racon round 2 polishing has already completed ***"
fi


#echo " *** generating overlap information after round 2 racon ***"

#if test ! -e results/${sample}/minimap2_nanopore_overlap3.sam.gz && test results/${sample}/racon_round2_output.fasta -nt results/${sample}/minimap2_nanopore_overlap3.sam.gz; then
#        echo Running minimap2 round 3
#        minimap2 -a -x ava-ont -t 32 results/${sample}/racon_round2_output.fasta ${nanopore_reads} > results/${sample}/minimap2_nanopore_overlap3.sam; gzip results/${sample}/minimap2_nanopore_overlap3.sam.gz
#		echo "*** Minimap2 overlap round 3 completed ***"
#else
#        echo "*** Minimap2 overlap round 3 has already completed ***"
#fi

#if test ! -e results/${sample}/racon_round3_output.fasta &&  test  results/${sample}/minimap2_nanopore_overlap3.sam.gz -nt results/${sample}/racon_round3_output.fasta; then
#        echo Running racon round 3
#        racon --threads 32 --include-unpolished ${nanopore_reads} results/${sample}/minimap2_nanopore_overlap3.sam.gz results/${sample}/racon_round2_output.fasta > results/${sample}/racon_round3_output.fasta
#		echo " **** Racon round 3 polishing completed ***"
#else
#        echo " **** Racon round 3 polishing has already completed ***"
#fi


## next step map illumina reads to Nanopore polished assembly from above steps

source bowtie2-2.2.9
source pilon-1.22
source samtools-1.9

echo "*** running bowtie and pilon ***"

raconfile=racon_round2_output.fasta
if test ! -e $raconfile; then
        bowtie2-build -f results/${sample}/$raconfile results/${sample}/$raconfile
	echo "*** bowtie2 build completed***"
fi


if [ "$illumina_reads_2" == "" && "$R3" == "" ]; then
	if test ! -e result/${sample}/illumina_aligned_to_assembly_sorted.bam; then
        	bowtie2 --threads 8 --qc-filter --no-unal -x results/${sample}/$raconfile -U $illumina_reads_1 $illumin_reads_2 $R3 $R4| samtools view -b  | samtools sort -o results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam; samtools index results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam
	else
		echo "*** bowtie has already completed***"
	fi

	if test ! -e results/${sample}/${sample}_pilon_polished_round1.fasta; then
		pilon --genome results/${sample}/$raconfile --unpaired results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam --outdir results/${sample} --output ${sample}_pilon_polished_round1.fasta --changes --threads 16
   	 else
		echo "*** Pilcon round 1 has already completed ***"
	fi

	# round 2 pilon polishing

	if test ! -e results/${sample}_pilon_polished_round2.fasta; then
		bowtie2 --threads 8 --qc-filter --no-unal -x results/${sample}_pilon_polished_round1.fasta -U $illumina_reads_1 $illumina_reads_2 $R3 $R4| samtools view -b  | samtools sort -o results/${sample}/illumina_aligned_to_assembly_sorted2.sorted.bam; samtools index results/${sample}/illumina_aligned_to_assembly_sorted2.sorted.bam
		echo "running second round of pilon"
		pilon --genome results/${sample}/${sample}_pilon_polished_round1.fasta --unpaired results/${sample}/illumina_aligned_to_assembly_sorted2.sorted.bam --outdir results/${sample} --output ${sample}_pilon_polished_round2.fasta -changes --threads 16
	fi

else
	if test ! -e results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam; then
		bowtie2 --threads 8 --qc-filter --no-unal -x results/${sample}/$raconfile  -1 $illumina_reads_1 -2 $illumina_reads_2 | samtools view -b  | samtools sort -o results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam; samtools index results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam;
	else
		echo "*** Bowtie has already completed ***"
	fi

	if test ! -e results/${sample}/${sample}_pilon_polished_round1.fasta; then
		pilon --genome results/${sample}/$raconfile --frags results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam --outdir results/${sample} --output ${sample}_pilon_polished_round1.fasta --changes --threads 16
	else
		echo "*** Pilcon round 1 has already completed ***"
    	fi
	# round 2 pilon polishing

	if test ! -e results/${sample}_pilon_polished_round2.fasta; then
		bowtie2 --threads 8 --qc-filter --no-unal -x results/${sample}/${sample}_pilon_polished_round1.fasta  -1 $illumina_reads_1 -2 $illumina_reads_2 | samtools view -b  | samtools sort -o results/${sample}/illumina_aligned_to_assembly_sorted2.sorted.bam; samtools index results/${sample}/illumina_aligned_to_assembly_sorted2.sorted.bam;
		echo "running second round of pilon"
		pilon --genome results/${sample}/results/${sample}/${sample}_pilon_polished_round1.fasta --frags results/${sample}/illumina_aligned_to_assembly_sorted2.sorted.bam --outdir results/${sample} --output ${sample}_pilon_polished_round2.fasta --changes --threads 16
	fi

fi
