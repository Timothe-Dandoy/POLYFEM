#!/bin/bash -l
#SBATCH -J Depth_Haplotype_EDNRB
#SBATCH -o ../../logs/depth_haplo_ednrb.log
#SBATCH -e ../../logs/depth_haplo_ednrb.log
#SBATCH -t 24:00:00
#SBATCH -p workq
#SBATCH --cpus-per-task 8
#SBATCH --mem=100G

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variables

#Need to give the name of the species you want to study => Ex: P_guadarramae, P_hispanicus, ...
species=$1
nickname=$(echo ${species} | cut -d'_' -f2 | cut -c1-3)
Name_Chr=NC_041330.1       # Chromosome Z
Region=30556153-30574830   # EDNRNB gene
diagnostic_region=30562000-30556900

#Inputs
file_pheno=../../../data/Pheno.xlsx
bam_list=../../../tmp/${species}/Files_pheno/${nickname}_bam_list_filtered.fofn

#Temporary
file_depth=../../../tmp/${species}/Files_genet/${nickname}_depth_edrnb.txt

#output
prefix_results=../../../results/${species}/Depth_EDNRB/${nickname}_plot_depth_ednrb_

#Secondary script
analyse_depth_R=../Secondary/analyse_depth.R

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Get depth EDNRB
echo "step 1"

module purge
module load bioinfo/samtools/1.23

> "$file_depth"

samtools depth -r ${Name_Chr}:${Region} -a -f "$bam_list" -o "$file_depth"

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2. Analyse depth 
echo "step 2"

module purge
module load statistics/R/4.5.0
module load compilers/gcc/12.2.0

Rscript $analyse_depth_R $file_depth $bam_list $file_pheno $prefix_results $diagnostic_region $species