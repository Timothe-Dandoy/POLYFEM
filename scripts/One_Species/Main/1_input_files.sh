#!/bin/bash -l
#SBATCH -J input_file
#SBATCH -o ../../logs/input_file.out
#SBATCH -e ../../logs/input_file.out
#SBATCH -p workq
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --mail-type=FAIL
#SBATCH --time=24:00:00

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variables

#Need to give the name of the species you want to study => Ex: P_guadarramae, P_hispanicus, ...
species=$1
nickname=$(echo ${species} | cut -d'_' -f2 | cut -c1-3)

#Inputs
file_pheno=../../../data/Pheno.xlsx
dir_tmp=../../../tmp/${species}/Files_pheno/${nickname}_ 
dir_sec=../Secondary/

#Outputs
bam_list=${dir_tmp}bam_list_filtered.fofn
bam_list_S=${dir_tmp}bam_list_filtered_S.fofn
bam_list_RB=${dir_tmp}bam_list_filtered_RB.fofn
individual_file=${dir_tmp}individual.txt
ind_pheno_file=${dir_tmp}individual_phenotype.txt
pop_file=${dir_tmp}population.file
pheno_file=${dir_tmp}phenotype.ybin

#Secondary scripts
build_input_R=${dir_sec}build_input_files.R

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. build input files in R
echo
echo "step 1"
echo

module purge
module load statistics/R/4.5.0
module load compilers/gcc/15.1.0

Rscript $build_input_R $file_pheno $species $bam_list $bam_list_S $bam_list_RB $individual_file $ind_pheno_file $pop_file $pheno_file

