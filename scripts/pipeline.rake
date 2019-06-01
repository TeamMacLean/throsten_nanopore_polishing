#!/usr/bin/env rake

ENV["sample"] ? @sample=ENV["sample"] : nil
ENV["nanopore_reads"] ? @nanopore_reads=ENV["nanopore_reads"] : nil
ENV["assembly"] ? @assembly=ENV["assembly"] : nil
ENV["R1"] ? @R1=ENV["R1"] : nil
ENV["R2"] ? @R2=ENV["R2"] : nil
ENV["R3"] ? @R3=ENV["R3"] : nil
ENV["R4"] ? @R4=ENV["R4"] : nil

directory "results/#{@sample}"

# minimap round 1
file "results/#{@sample}/minimap2_nanopore_overlap1.sam.gz" =>  ["results/#{@sample}", "#{@assembly}", "#{@nanopore_reads}"] do
    puts "Running minimap round 1"
    sh "source minimap2-2.13; minimap2 -a -x ava-ont -t 32 #{@assembly} #{@nanopore_reads} | gzip > results/${sample}/minimap2_nanopore_overlap1.sam.gz"
end

# racon round 1
file "results/#{@sample}/racon_round1_output.fasta" => ["results/#{@sample}/minimap2_nanopore_overlap1.sam.gz"] do
    puts "Running racon roudn 1"
    sh "source racon-1.3.2;  racon --threads 32 --include-unpolished  #{@nanopore_reads} results/#{@sample}/minimap2_nanopore_overlap1.sam.gz #{@assembly} > results/#{@sample}/racon_round1_output.fasta"
end
# minimap round2

file "results/#{@sample}/minimap2_nanopore_overlap2.sam.gz" => ["results/#{@sample}/racon_round1_output.fasta"] do
    puts "Running minimap round 2"
    sh "source minimap2-2.13;  minimap2 -a -x ava-ont -t 32 results/#{@sample}/racon_round1_output.fasta #{@nanopore_reads}  | gzip > results/#{@sample}/minimap2_nanopore_overlap2.sam.gz"
end

# racond round 2
file "results/#{@sample}/racon_round2_output.fasta" => [ "results/#{@sample}/minimap2_nanopore_overlap2.sam.gz", "results/#{@sample}/racon_round1_output.fasta"] do
    puts "Running racon round 2"
    sh "source racon-1.3.2;  racon --threads 32 --include-unpolished #{@nanopore_reads} results/#{@sample}/minimap2_nanopore_overlap2.sam.gz results/#{@sample}/racon_round1_output.fasta > results/#{@sample}/racon_round2_output.fasta"
end

#pilon round 1 for FR13, BR32, US71
file "results/#{@sample}/racon_round2_output.fasta.1.bt2" => "results/#{@sample}/racon_round2_output.fasta" do
    puts "Building bowtie index for pilon round 1"
    sh "./scripts/build_bowtie_index.sh results/#{@sample}/racon_round2_output.fasta"
end
file "results/#{@sample}/illumina_aligned_to_assembly_sorted.sorted.bam" => ["results/#{@sample}/racon_round2_output.fasta.1.bt2"] do
    puts "running bowtie2 for pilon round 1"
    sh "source bowtie2-2.2.9; source samtools-1.9; bowtie2 --threads 2 --qc-filter --no-unal -x results/#{@sample}/racon_round2_output.fasta -U #{@R1} -S results/#{@sample}/illumina_aligned_to_assembly1.sam && samtools view -b results/#{@sample}/illumina_aligned_to_assembly1.sam | samtools sort -o results/#{@sample}/illumina_aligned_to_assembly_sorted.sorted.bam && samtools index results/#{@sample}/illumina_aligned_to_assembly_sorted.sorted.bam && rm results/#{@sample}/illumina_aligned_to_assembly1.sam"
end

file "results/#{@sample}/#{@sample}_pilon_polished_round1.fasta" => "results/#{@sample}/illumina_aligned_to_assembly_sorted.sorted.bam" do
    puts "running pilon round 1"
    sh "source pilon-1.22;  pilon --genome results/#{@sample}/racon_round2_output.fasta --frags results/#{@sample}/illumina_aligned_to_assembly_sorted.sorted.bam --outdir results/#{@sample}/ --output #{@sample}_pilon_polished_round1 --changes --threads 16 "
end

# pilon round 2 polishing
file "results/#{@sample}/#{@sample}_pilon_polished_round1.fasta.1.bt2" => "results/#{@sample}/#{@sample}_pilon_polished_round1.fasta" do
    puts "Building index fro pilon poslished round1 fasta "
    sh "./scripts/build_bowtie_index.sh results/#{@sample}/#{@sample}_pilon_polished_round1.fasta"
end

file "results/#{@sample}/illumina_aligned_to_assembly_sorted2.sorted.bam" => [ "results/#{@sample}/#{@sample}_pilon_polished_round1.fasta.1.bt2" , "#{@R1}" ] do
    puts "Aligning reads to pilon polished round1 fasta"
    sh "source bowtie2-2.2.9; source samtools-1.9;  bowtie2 --threads 8 --qc-filter --no-unal -x results/#{@sample}/#{@sample}_pilon_polished_round1.fasta -U #{@R1} -S results/#{@sample}/illumina_aligned_to_assembly2.sam && samtools view -b results/#{@sample}/illumina_aligned_to_assembly2.sam | samtools sort -o results/#{@sample}/illumina_aligned_to_assembly_sorted2.sorted.bam && samtools index results/#{@sample}/illumina_aligned_to_assembly_sorted2.sorted.bam && rm results/#{@sample}/illumina_aligned_to_assembly2.sam"
end

file "results/#{@sample}/#{@sample}_pilon_polished_round2.fasta" => [ "results/#{@sample}/illumina_aligned_to_assembly_sorted2.sorted.bam"] do
    sh "source pilon-1.22;  pilon --genome results/#{@sample}/#{@sample}_pilon_polished_round1.fasta --frags results/#{@sample}/illumina_aligned_to_assembly_sorted2.sorted.bam --outdir results/#{@sample} --output #{@sample}_pilon_polished_round2 --changes --threads 16"
end

#pilon round 1 for CD156
file "results/#{@sample}/#{@sample}_illumina_aligned_to_assembly_sorted.sorted.bam" => ["results/#{@sample}/racon_round2_output.fasta.1.bt2", "#{@R1}", "#{@R2}", "#{@R3}", "#{@R4}"] do
    puts "running bowtie2 for pilon round 1 for CD156 "
    sh "source bowtie2-2.2.9; source samtools-1.9;  bowtie2 --threads 2 --qc-filter --no-unal -x results/#{@sample}/racon_round2_output.fasta -U #{@R1} -U #{@R2} -U #{@R3} -U #{@R4} -S results/#{@sample}/illumina_aligned_to_assembly1.sam &&  samtools view -b results/#{@sample}/illumina_aligned_to_assembly1.sam | samtools sort -o results/#{@sample}/#{@sample}_illumina_aligned_to_assembly_sorted.sorted.bam && ln -s results/#{@sample}/#{@sample}_illumina_aligned_to_assembly_sorted.sorted.bam results/#{@sample}/illumina_aligned_to_assembly_sorted.sorted.bam && samtools index results/#{@sample}/illumina_aligned_to_assembly_sorted.sorted.bam && rm results/#{@sample}/illumina_aligned_to_assembly1.sam"
end

file "results/#{@sample}/#{@sample}_illumina_aligned_to_assembly_sorted2.sorted.bam" => [ "results/#{@sample}/#{@sample}_pilon_polished_round1.fasta.1.bt2" , "#{@R1}", "#{@R2}", "#{@R3}", "#{@R4}" ] do
    puts "Aligning reads to pilon polished round1 fasta for CD156"
    sh "source bowtie2-2.2.9; source samtools-1.9; bowtie2 --threads 8 --qc-filter --no-unal -x results/#{@sample}/#{@sample}_pilon_polished_round1.fasta -U #{@R1} -U #{@R2} -U #{@R3} -U #{@R4} -S results/#{@sample}/illumina_aligned_to_assembly2.sam && samtools view -b results/#{@sample}/illumina_aligned_to_assembly2.sam | samtools sort -o results/#{@sample}/#{@sample}_illumina_aligned_to_assembly_sorted2.sorted.bam && ln -s results/#{@sample}/#{@sample}_illumina_aligned_to_assembly_sorted2.sorted.bam results/#{@sample}/illumina_aligned_to_assembly_sorted2.sorted.bam && samtools index results/#{@sample}/illumina_aligned_to_assembly_sorted2.sorted.bam && rm results/#{@sample}/illumina_aligned_to_assembly2.sam"
end

task :CD156 => ["results/#{@sample}/racon_round2_output.fasta.1.bt2", "results/#{@sample}/#{@sample}_illumina_aligned_to_assembly_sorted.sorted.bam", "results/#{@sample}/#{@sample}_pilon_polished_round1.fasta", "results/#{@sample}/#{@sample}_illumina_aligned_to_assembly_sorted2.sorted.bam","results/#{@sample}/#{@sample}_pilon_polished_round2.fasta"  ] do
    puts "running sample CD156"
end
task :Other => ["results/#{@sample}/racon_round2_output.fasta.1.bt2", "results/#{@sample}/illumina_aligned_to_assembly_sorted.sorted.bam", "results/#{@sample}/#{@sample}_pilon_polished_round1.fasta", "results/#{@sample}/illumina_aligned_to_assembly_sorted2.sorted.bam", "results/#{@sample}/#{@sample}_pilon_polished_round2.fasta" ] do
    puts "running sample #{@sample}"
end


task :default do
    puts "#{@sample} #{@nanopore_reads} #{@assembly}"
    if "#{@sample}" == "CD156"
        Rake::Task[:CD156].invoke
    else
        puts "sample #{@sample}"
        Rake::Task[:Other].invoke
    end

end
