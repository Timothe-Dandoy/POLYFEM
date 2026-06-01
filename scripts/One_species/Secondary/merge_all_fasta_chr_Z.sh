#!/bin/bash -l
#SBATCH -J Merge_all_fasta_chr_Z
#SBATCH -o ../../logs/merge_fasta_z.out
#SBATCH -e ../../logs/merge_fasta_z.out
#SBATCH -t 24:00:00
#SBATCH -p workq
#SBATCH --cpus-per-task 8
#SBATCH --mem=100G

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variables

species=$1
nickname=$(echo ${species} | cut -d'_' -f2 | cut -c1-3)
Name_Chr=NC_041330.1

prefix_fastas_chr_aligned=../../../tmp/${species}/Files_genet/fasta/${nickname}_${Name_Chr}_allind_
sufix_fastas_chr_aligned=_aligned.fa 

fasta_chr_merged=../../../tmp/${species}/Files_genet/fasta/${nickname}_${Name_Chr}_allind_merged.fa  

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Merge fasta Chromosome Z
echo ""
echo "step 1"
echo ""

module load bioinfo/SeqKit/2.12.0

seqkit concat $(ls ${prefix_fastas_chr_aligned}*${sufix_fastas_chr_aligned}* | sort -V) > $fasta_chr_merged
