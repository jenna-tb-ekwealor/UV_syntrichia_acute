#!/bin/bash
# Job name:
#SBATCH --job-name=acute_paired-end_tophat
#
# Account:
#SBATCH --account=fc_phylodiv
#
# Partition:
#SBATCH --partition=savio
#
# Quality of Service:
#SBATCH --qos=savio_normal
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
## Command(s) to run:

export PATH=$PATH:/global/home/groups/fc_phylodiv/modules/
export PATH=$PATH:/global/home/groups/fc_phylodiv/modfiles/

# run from "clean" directory

module load bowtie2

for sample in */; do
 	samplebname=$(basename "$sample")
	mkdir tophat_out_"$samplebname"
	/global/home/groups/fc_phylodiv/modules/tophat/tophat-2.1.1.Linux_x86_64/tophat -r 275 -o tophat_out_"$samplebname" --no-coverage-search --transcriptome-index /global/scratch/jbaughman/caninervis_genome_assembly/bowtie_index/s_caninervis/transcriptome_index /global/scratch/jbaughman/caninervis_genome_assembly/bowtie_index/s_caninervis "$sample"*_1_clean.fastq "$sample"*_2_clean.fastq 
 done