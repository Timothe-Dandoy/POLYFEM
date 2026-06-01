#!/bin/bash -l
#SBATCH -J input_file_all_species
#SBATCH -o ../../logs/input_file_allspe.out
#SBATCH -e ../../logs/input_file_allspe.out
#SBATCH -p workq
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --mail-type=FAIL
#SBATCH --time=24:00:00

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variables

file_pheno=../../../data/Pheno.xlsx
dir_tmp=../../../tmp/All_species/Files_pheno/
dir_sec=../Secondary/

bam_list=${dir_tmp}bam_allspecies_female_filtered.fofn

build_input_R=${dir_sec}build_input_files.R

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. build input files in R
echo
echo "step 1"
echo

module purge
module load statistics/R/4.5.0
module load compilers/gcc/15.1.0

Rscript $build_input_R $file_pheno $bam_list

