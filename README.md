## Introduction

We are re-attempting to polish  nanopore data assembly using minimap2, racon and pilon with nanopore long reads and illumina short reads. We will do two rounds of minimap2 and racon step as it is suggested that the assembly polishing improves while repeating the minimap2 and racon step. Finally, the nanopore polished assembly with minimap2 and racon is again repolished with two rounds of pilon using nextgen illumina short reads.

The scripts are named for each sample. For the second attempt to polish the assemblies, the scripts are named as version2. All scripts have absolute paths to filenames. Any change in input files means the filenames will need to be manually added in the script.
1) run_fr13_version2.sh - runs sample FR13 assembly polishing
2) run_cd156_version2.sh - runs sample CD156 assembly polishing
3) run_br32_version2.sh - runs sample BR32 assembly polishing
4) run_us71_version2.sh - runs sample US71 assembly polishing


