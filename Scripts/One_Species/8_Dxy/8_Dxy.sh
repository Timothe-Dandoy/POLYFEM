#!/bin/bash -l
#SBATCH -J Genetic_distance
#SBATCH -o ../../logs/genetic_distance.log
#SBATCH -e ../../logs/genetic_distance.log
#SBATCH -t 24:00:00
#SBATCH -p workq
#SBATCH --cpus-per-task 8
#SBATCH --mem=500G

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variables

#Data
# Need to give the name of the species you want to study => Ex: P_guadarramae, P_hispanicus, ...
species=$1
nickname=$(echo ${species} | cut -d'_' -f2 | cut -c1-3) 
Name_Chr=NC_041330.1   # chromosome Z
diagnostic_region=30556153-30574830 # whole EDNRB gene

#Inputs
ind_pheno=../../../tmp/${species}/Files_pheno/${nickname}_individual_phenotype.txt
fasta_allind_ednrb_aligned=../../../tmp/${species}/Files_genet/${nickname}_allind_ednrb_aligned.fa
fasta_chr_merged=../../../tmp/${species}/Files_genet/fasta/${nickname}_${Name_Chr}_allind_merged.fa  

#Outputs
prefix_output_ednrb=../../../results/${species}/Dxy/${nickname}_ednrb_dxy_all_couple
prefix_output_chrZ=../../../results/${species}/Dxy/${nickname}_${Name_Chr}_chrz_dxy_pheno_ 

#Secondary scripts
distrib_dxy_ednrb_R=../Secondary/distribution_dxy_ednrb_allpair.R
compute_dxy_chrZ_R=../Secondary/compute_dxy_chrZ_pheno.R

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Compute genetic distance along the EDNRB sequence
echo ""
echo "step 1"
echo ""

module purge
module load statistics/R/4.5.0
module load compilers/gcc/15.1.0

Rscript $distrib_dxy_ednrb_R $fasta_allind_ednrb_aligned $ind_pheno $species $prefix_output_ednrb 


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2. Compute genetic distance along the chromosome Z
echo ""
echo "step 2"
echo ""

module purge  
module load statistics/R/4.5.0
module load compilers/gcc/15.1.0

Rscript $compute_dxy_chrZ_R $fasta_chr_merged $ind_pheno $diagnostic_region $species $prefix_output_chrZ 

