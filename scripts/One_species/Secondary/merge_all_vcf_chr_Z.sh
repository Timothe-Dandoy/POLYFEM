#!/bin/bash -l
#SBATCH -J merge_VCF_chr_Z
#SBATCH -o ../../logs/merge_vcf_z.out
#SBATCH -e ../../logs/merge_vcf_z.out
#SBATCH -p workq
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --mail-type=FAIL

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variables

species=$1
nickname=$(echo ${species} | cut -d'_' -f2 | cut -c1-3)
chr=NC_041330.1
dir_tmp=../../../tmp/${species}/Files_genet/vcf/

START_IDX=1   
END_IDX=509     

merged_tmp=${dir_tmp}${nickname}_${chr}_vcf_unsorted.vcf.gz
merged_sorted=${dir_tmp}${nickname}_${chr}_vcf_merged.vcf.gz

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Generate the list of VCFs
echo
echo "step 1"
echo

vcf_list=$(ls ${dir_tmp}*.vcf | sort -V | sed -n "${START_IDX},${END_IDX}p")

# Ensure there are enough files
if [[ -z ${vcf_list} ]]; then
  echo "No VCF files found in the specified range: $START_IDX to $END_IDX"
  exit 1
fi

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2. Conctenate, sort and index
echo
echo "step 2"
echo

module purge
module load bioinfo/Bcftools

bcftools concat -Oz -o ${merged_tmp} ${vcf_list}
bcftools sort -Oz ${merged_tmp} -o ${merged_sorted}
bcftools index ${merged_sorted}

rm ${merged_tmp}

