#!/bin/bash

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variables

# Need to give the name of the species you want to study => Ex: P_guadarramae, P_hispanicus, ...
species=$1

#Secondary scripts
build_vcf_bash=../Secondary/build_region_vcf_chr_Z.sh
merge_vcf_bash=../Secondary/merge_all_vcf_chr_Z.sh

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Build VCF for all region of chr Z
echo
echo "step 1"
echo

ARRAY_JOB_ID=$(sbatch --parsable ${build_vcf_bash} ${species})

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2. Merge all vcf
echo
echo "step 2"
echo

sbatch --dependency=afterok:${ARRAY_JOB_ID} ${merge_vcf_bash} ${species}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3. Remove .out
echo
echo "step 3"
echo

rm *slurm*