#!/bin/bash -l
#SBATCH -J chr_Z_LD
#SBATCH -o ../../logs/ld_z.log
#SBATCH -e ../../logs/ld_z.log
#SBATCH -t 01:00:00
#SBATCH -p workq
#SBATCH --cpus-per-task 8
#SBATCH --mem=10G

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variables

#Data
# Need to give the name of the species you want to study => Ex: P_guadarramae, P_hispanicus, ...
species=$1
nickname=$(echo ${species} | cut -d'_' -f2 | cut -c1-3)
CHR=NC_041330.1  # Chromosome Z 
START_IDX=1    # indice +1
END_IDX=509     # indice +1
window=200
window_mean=5000
diagnostic_region=30556153-30574830  # whole EDNRB gene

#Inputs
file_vcf=../../../tmp/${species}/Files_genet/vcf/${nickname}_${CHR}_vcf_merged.vcf.gz

#Temporary
prefix_plink=../../../tmp/${species}/Files_genet/${nickname}_${CHR}_plink
prefix_plink_window=../../../tmp/${species}/Files_genet/${nickname}_${CHR}_plink_${window}bp
file_ld=../../../tmp/${species}/Files_genet/${nickname}_${CHR}_plink_${window}bp.ld
file_bim=../../../tmp/${species}/Files_genet/${nickname}_${CHR}_plink_${window}bp.bim

#Output
prefix_plot=../../../results/${species}/LD/${nickname}_${CHR}_plink_${window}bp_LD

# Secondary scripts
analyse_ld=../Secondary/analyse_ld_plink.R

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. build plink format
echo ""
echo "step 1"
echo ""

module purge
module load bioinfo/PLINK/1.90b7

# plink --vcf $file_vcf \
#       --allow-extra-chr \
#       --double-id \
#       --out $prefix_plink
      
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2. Filtering for 1 SNP each 200 bp 
echo
echo "step 2"
echo

module purge
module load bioinfo/PLINK/1.90b7

# plink --bfile $prefix_plink \
#       --allow-extra-chr \
#       --bp-space $window \
#       --make-bed \
#       --out $prefix_plink_window

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3. Compute LD between SNP 
echo
echo "step 3"
echo

module purge
module load bioinfo/PLINK/1.90b7

# plink \
#      --bfile $prefix_plink_window \
#      --r2 \
#      --allow-extra-chr \
#      --ld-window 2 \
#      --ld-window-kb 1000\
#      --ld-window-r2 0 \
#      --out $prefix_plink_window

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 4. Plot LD along the Z chromosome
echo
echo "step 4"
echo

module purge
module load statistics/R/4.5.0
module load compilers/gcc/15.1.0

Rscript $analyse_ld $file_ld $file_bim $diagnostic_region $window_mean $species $prefix_plot 