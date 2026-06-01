#!/bin/bash -l
#SBATCH -J GWAS
#SBATCH -o ../../logs/GWAS.out
#SBATCH -e ../../logs/GWAS.out
#SBATCH -p workq
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G
#SBATCH --mail-type=FAIL
#SBATCH --time=72:00:00

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variables

species=$1
nickname=$(echo ${species} | cut -d'_' -f2 | cut -c1-3)

File_Pheno=../../../tmp/${species}/Files_pheno/${nickname}_phenotype.ybin
File_Pop=../../../tmp/${species}/Files_pheno/${nickname}_population.txt
File_BamList=../../../tmp/${species}/Files_pheno/${nickname}_bam_list_filtered.fofn
prefix_out=../../../tmp/${species}/Files_genet/${nickname}_out_gwas
out_GWAS=../../../tmp/${species}/Files_genet/${nickname}_out_gwas.lrt0.gz

prefix_res=../../../results/${species}/GWAS/${nickname}_gwas_

Analyse_R=../Secondary/Analyse_GWAS.R

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. GWAS Computation
echo
echo "step 1"
echo

module load bioinfo/ANGSD/0.940

angsd -yQuant $File_Pheno \
      -doAsso 2 \
      -cov $File_Pop \
      -GL 1 \
      -doPost 1 \
      -out $prefix_out \
      -doMajorMinor 1 \
      -doMaf 1 \
      -SNP_pval 1e-6 \
      -Pvalue 1 \
      -bam $File_BamList \
      -P 8 

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2. GWAS Analysis
echo
echo "step 2"
echo

module load statistics/R/4.5.0
module load compilers/gcc/15.1.0

Rscript $Analyse_R $out_GWAS $species $prefix_res



