#!/bin/bash -l
#SBATCH -J build_VCF_chr_Z
#SBATCH -o ../../logs/vcf/build_vcf_z_%a.out
#SBATCH -e ../../logs/vcf/build_vcf_z_%a.out
#SBATCH -p unlimitq
#SBATCH --array=0-508
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --mail-type=FAIL

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variables

#Data
species=$1
nickname=$(echo ${species} | cut -d'_' -f2 | cut -c1-3)
CHR=NC_041330.1        # Chromosome Z

#Inputs
REGION_FILE=../../../data/P_muralis/${CHR}.regions.txt
bam_list=../../../tmp/${species}/Files_pheno/${nickname}_bam_list_filtered.fofn
ref=../../../../ref/GCF_004329235.1_PodMur_1.0_genomic.fna
dir_tmp=../../../tmp/${species}/Files_genet/vcf/

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Get the sub region for this SLURM job
echo
echo "step 1"
echo

REGION=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ${REGION_FILE})

if [ -z ${REGION} ]; then
    echo "Error: REGION is empty for task ID $SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2. Build the VCF for the sub region
echo
echo "step 2"
echo

module purge
module load bioinfo/freebayes/1.3.6

freebayes -r ${REGION} -f ${ref} -L ${bam_list} \
  --min-coverage 100 --limit-coverage 100 -E -1 --max-complex-gap -1 \
  --haplotype-length -1 --hwe-priors-off -m 20 -q 20 -p 1 -n 4 \
  > "${dir_tmp}region_${SLURM_ARRAY_TASK_ID}.vcf"
