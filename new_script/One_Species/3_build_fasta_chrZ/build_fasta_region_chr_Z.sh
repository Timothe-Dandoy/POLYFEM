#!/bin/bash -l
#SBATCH -J Build_fasta_chr_Z
#SBATCH -o ../../logs/fasta/build_fasta_z_%a.out
#SBATCH -e ../../logs/fasta/build_fasta_z_%a.out
#SBATCH -t 72:00:00
#SBATCH --array=2-102%20
#SBATCH -p workq
#SBATCH --cpus-per-task 8
#SBATCH --mem=50G

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variables

#Data
species=$1
nickname=$(echo ${species} | cut -d'_' -f2 | cut -c1-3)
Name_Chr=NC_041330.1     # Chromosome Z
Reg_studied="$((1+SLURM_ARRAY_TASK_ID*500000))"-"$(((SLURM_ARRAY_TASK_ID+1)*500000))"

if [ $SLURM_ARRAY_TASK_ID -eq 102 ]
    then Reg_studied=50500001-50871175
fi

#Inputs
Bam_list=../../../tmp/${species}/Files_pheno/${nickname}_bam_list_filtered.fofn
fasta_ref_chr=../../../data/P_muralis/ref_chrZ.fa

#Temporary
fasta_region=../../../tmp/${species}/Files_genet/fasta/${nickname}_${Name_Chr}_allind_${SLURM_ARRAY_TASK_ID}.fa

#Output
fasta_region_aligned=../../../tmp/${species}/Files_genet/fasta/${nickname}_${Name_Chr}_allind_${SLURM_ARRAY_TASK_ID}_aligned.fa 

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Get fasta of the region
echo 
echo "step 1"
echo 

module purge
module load bioinfo/samtools/1.23
 
> $fasta_region

while read line
do
    File=$line
    Individual=$(basename $File .bam)
    samtools consensus -r ${Name_Chr}:${Reg_studied} -T $fasta_ref_chr -a -f fasta $File | \
    sed "s/^>.*/>$Individual/" >> $fasta_region
done < $Bam_list

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2. Aligne Sequences of every individual
echo
echo "step 2"
echo

module purge
module load bioinfo/MAFFT/7.505

mafft --auto --thread $SLURM_CPUS_PER_TASK ${fasta_region} > ${fasta_region_aligned} 

