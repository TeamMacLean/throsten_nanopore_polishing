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

if test ! -e results/${sample}/racon_round1_output.fasta ||  test results/${sample}/minimap2_nanopore_overlap1.sam.gz -nt results/${sample}/racon_round1_output.fasta || test ! -s  results/${sample}/racon_round1_output.fasta; then
        echo "Running racon round 1"
        racon --threads 32 --include-unpolished  ${nanopore_reads} results/${sample}/minimap2_nanopore_overlap1.sam.gz ${nanopore_assembly} > results/${sample}/racon_round1_output.fasta
	echo "*** Racon round 1 polishing completed ***"
else
        echo "*** Racon round 1 polishing has already completed ***"
fi

echo "*** generating overlap information after round 1 racon ***"

if test ! -e results/${sample}/minimap2_nanopore_overlap2.sam.gz || test results/${sample}/racon_round1_output.fasta -nt results/${sample}/minimap2_nanopore_overlap2.sam.gz; then
        echo Running minimap2 round 2
        minimap2 -a -x ava-ont -t 32 results/${sample}/racon_round1_output.fasta ${nanopore_reads}  | gzip > results/${sample}/minimap2_nanopore_overlap2.sam.gz
	echo "*** Minimap2 overlap round 2 completed ***"
else
        echo "*** Minimap2 overlap round 2 has already completed ***"
fi

if test ! -e results/${sample}/racon_round2_output.fasta ||  test results/${sample}/minimap2_nanopore_overlap2.sam.gz -nt results/${sample}/racon_round2_output.fasta || test ! -s results/${sample}/racon_round2_output.fasta; then
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

function bowtiebuild(){
    if test ! -e $1.1.bt2; then
        bowtie2-build -f $1 $1
    	echo "*** bowtie2 build completed***"
    fi
}
echo "*** running bowtie and pilon ***"

raconfile=racon_round2_output.fasta



if [[ "$sample" == "FR13" || "$sample" == "US71" || "$sample" == "BR32" ]]; then
    echo "Running bowtie pilon for $sample"
    if test ! -e result/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam; then
        echo "bowtie2 --threads 8 --qc-filter --no-unal -x results/${sample}/$raconfile -U $illumina_reads_1 -S results/${sample}/illumina_aligned_to_assembly1.sam"
        cmd="bowtie2 --threads 8 --qc-filter --no-unal -x results/${sample}/$raconfile -U $illumina_reads_1 -S results/${sample}/illumina_aligned_to_assembly1.sam ; samtools view -b results/${sample}/illumina_aligned_to_assembly1.sam | samtools sort -o results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam; samtools index results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam"
        $cmd
    else
        echo "*** bowtie has already completed***"
    fi

    if test ! -e results/${sample}/${sample}_pilon_polished_round1.fasta; then
        pilon results/${sample}/$raconfile results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam results/${sample} ${sample}_pilon_polished_round1

   	 else
		echo "*** Pilcon round 1 has already completed ***"
	fi

	# round 2 pilon polishing

	if test ! -e results/${sample}_pilon_polished_round2.fasta; then
        bowtiebuild results/${sample}/${sample}_pilon_polished_round1.fasta
		bowtie2 --threads 8 --qc-filter --no-unal -x results/${sample}/${sample}_pilon_polished_round1.fasta -U $illumina_reads_1 -S results/${sample}/illumina_aligned_to_assembly2.sam; samtools view -b results/${sample}/illumina_aligned_to_assembly2.sam | samtools sort -o results/${sample}/illumina_aligned_to_assembly_sorted2.sorted.bam; samtools index results/${sample}/illumina_aligned_to_assembly_sorted2.sorted.bam
		echo "running second round of pilon"
		run-pilon results/${sample}/${sample}_pilon_polished_round1.fasta results/${sample}/illumina_aligned_to_assembly_sorted2.sorted.bam results/${sample} ${sample}_pilon_polished_round2
	fi

elif [[ "$sample" == "CD156" ]]; then
	echo "running bowtie pilon for $sample"
	if [ ! -e result/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam ]; then
		echo "bowtie2 --threads 8 --qc-filter --no-unal -x results/${sample}/$raconfile -U $illumina_reads_1 -U $illumina_reads_2 -U $R3 -U $R4 -S results/${sample}/illumina_aligned_to_assembly1.sam"
        cmd="bowtie2 --threads 8 --qc-filter --no-unal -x results/${sample}/$raconfile -U $illumina_reads_1 -U $illumina_reads_2 -U $R3 -U $R4 -S results/${sample}/illumina_aligned_to_assembly1.sam ; samtools view -b results/${sample}/illumina_aligned_to_assembly1.sam | samtools sort -o results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam; samtools index results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam"
        $cmd
	else
		echo "*** bowtie has already completed***"
	fi

	if test ! -e results/${sample}/${sample}_pilon_polished_round1.fasta; then
        echo "Running pilon round 1"
		run-pilon results/${sample}/$raconfile results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam results/${sample} ${sample}_pilon_polished_round1

   	 else
		echo "*** Pilcon round 1 has already completed ***"
	fi

	# round 2 pilon polishing

	if test ! -e results/${sample}_pilon_polished_round2.fasta; then
        echo running bowtie round 2
        bowtiebuild results/${sample}/${sample}_pilon_polished_round1.fasta
		bowtie2 --threads 8 --qc-filter --no-unal -x results/${sample}/${sample}_pilon_polished_round1.fasta -U $illumina_reads_1 -U $illumina_reads_2 -U $R3 -U $R4 -S results/${sample}/illumina_aligned_to_assembly2.sam; samtools view -b results/${sample}/illumina_aligned_to_assembly2.sam | samtools sort -o results/${sample}/illumina_aligned_to_assembly_sorted2.sorted.bam; samtools index results/${sample}/illumina_aligned_to_assembly_sorted2.sorted.bam
		echo "running second round of pilon"
		run-pilon results/${sample}/${sample}_pilon_polished_round1.fasta results/${sample}/illumina_aligned_to_assembly_sorted2.sorted.bam results/${sample} ${sample}_pilon_polished_round2
	fi

else
	echo "running bowtie pilon with paired data"
	if test ! -e results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam; then
		echo bowtie2 --threads 8 --qc-filter --no-unal -x results/${sample}/$raconfile  -1 $illumina_reads_1 -2 $illumina_reads_2 -S results/${sample}/illumina_aligned_to_assembly1.sam
		bowtie2 --threads 8 --qc-filter --no-unal -x results/${sample}/$raconfile  -1 $illumina_reads_1 -2 $illumina_reads_2 -S results/${sample}/illumina_aligned_to_assembly1.sam ; samtools view -b  results/${sample}/illumina_aligned_to_assembly1.sam | samtools sort -o results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam; samtools index results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam;
	else
		echo "*** Bowtie has already completed ***"
	fi

	if test ! -e results/${sample}/${sample}_pilon_polished_round1.fasta; then
		run-pilon  results/${sample}/$raconfile  results/${sample}/illumina_aligned_to_assembly_sorted.sorted.bam  results/${sample}  ${sample}_pilon_polished_round1
	else
		echo "*** Pilcon round 1 has already completed ***"
    fi
	# round 2 pilon polishing

	if test ! -e results/${sample}_pilon_polished_round2.fasta; then
        bowtiebuild results/${sample}/${sample}_pilon_polished_round1.fasta
		bowtie2 --threads 8 --qc-filter --no-unal -x results/${sample}/${sample}_pilon_polished_round1.fasta  -1 $illumina_reads_1 -2 $illumina_reads_2 -S results/${sample}/illumina_aligned_to_assembly2.sam; samtools view -b  results/${sample}/illumina_aligned_to_assembly2.sam | samtools sort -o results/${sample}/illumina_aligned_to_assembly_sorted2.sorted.bam; samtools index results/${sample}/illumina_aligned_to_assembly_sorted2.sorted.bam;
		echo "running second round of pilon"
		run-pilon results/${sample}/results/${sample}/${sample}_pilon_polished_round1.fasta  results/${sample}/illumina_aligned_to_assembly_sorted2.sorted.bam  results/${sample}  ${sample}_pilon_polished_round2
	fi

fi

function run-pilon(){
    if [ ! -e ${$3}/${4}.fasta ]; then
        pilon --genome $1 --frags $2 --outdir $3 --output $4 --changes --threads 16
    fi
}
