#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variables

# Need to give the name of the species you want to study => Ex: P_guadarramae, P_hispanicus, ...
species=$1

#Secondary files
build_fasta_bash=../Secondary/build_fasta_region_chr_Z.sh
merge_fasta_bash=../Secondary/merge_all_fasta_chr_Z.sh

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Build fasta for all region of chr Z
echo
echo "step 1"
echo

ARRAY_JOB_ID=$(sbatch --parsable ${build_fasta_bash} ${species})

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2. Merge all fasta
echo
echo "step 2"
echo

sbatch --dependency=afterok:${ARRAY_JOB_ID} ${merge_fasta_bash} ${species}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3. Remove .out
echo
echo "step 3"
echo
 
rm *slurm*