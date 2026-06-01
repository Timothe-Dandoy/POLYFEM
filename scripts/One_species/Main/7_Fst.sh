#!/bin/bash -l
#SBATCH -J angsd_fst_chr_Z
#SBATCH -o ../../logs/angsd_fst_chrz.out
#SBATCH -e ../../logs/angsd_fst_chrz.out
#SBATCH -p workq
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --mail-type=FAIL
#SBATCH --time=24:00:00

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variables

species=$1
nickname=$(echo ${species} | cut -d'_' -f2 | cut -c1-3)
CHR=NC_041330.1
window_size=5000
window_step=5000
diagnostic_region=30556153-30574830 # whole EDNRB gene

fasta_ref=../../../../ref/GCF_004329235.1_PodMur_1.0_genomic.fna 
bam_list_S=../../../tmp/${species}/Files_pheno/${nickname}_bam_list_filtered_S.fofn
bam_list_RB=../../../tmp/${species}/Files_pheno/${nickname}_bam_list_filtered_RB.fofn

prefix_saf_S=../../../tmp/${species}/Files_genet/${nickname}_saf_S
prefix_saf_RB=../../../tmp/${species}/Files_genet/${nickname}_saf_RB
twodsfs_prior=../../../tmp/${species}/Files_genet/${nickname}_prior_S_RB.ml
fst_out=../../../tmp/${species}/Files_genet/${nickname}_fst_S_RB.out
sliding_window_out=../../../tmp/${species}/Files_genet/${nickname}_fst_S_RB_sliding_window.out
prefix_plot=../../../results/${species}/Fst/${nickname}_fst_S_RB_sliding_window

tmp_file=../../../tmp/${species}/Files_genet/${nickname}_tmp_fst.out

plot_fst_R=../Secondary/analyse_fst_angsd.R

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Calculate FST
echo
echo "step 1"
echo

module purge
module load bioinfo/ANGSD/0.940

# # #first calculate per pop saf for each population
# angsd -bam ${bam_list_S}  -anc ${fasta_ref} -out ${prefix_saf_S} -dosaf 1 -gl 1 -r ${CHR}
# angsd -bam ${bam_list_RB}  -anc ${fasta_ref} -out ${prefix_saf_RB} -dosaf 1 -gl 1 -r ${CHR}

# #calculate the 2dsfs prior
# realSFS ${prefix_saf_S}.saf.idx ${prefix_saf_RB}.saf.idx > ${twodsfs_prior}

# #prepare the fst for easy window analysis
# realSFS fst index ${prefix_saf_S}.saf.idx ${prefix_saf_RB}.saf.idx -sfs ${twodsfs_prior} -fstout ${fst_out}

# #get the global estimate
# realSFS fst stats ${fst_out}.fst.idx

# #get the sliding windows estimate
# realSFS fst stats2 ${fst_out}.fst.idx -win ${window_size} -step ${window_step} > ${sliding_window_out}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2. Change the header
echo
echo "step 2"
echo

# sed '1c\region\tchr\tmidPos\tNsites\tfst' ${sliding_window_out} > ${tmp_file}

# mv ${tmp_file} ${sliding_window_out}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3. Plot Fst
echo
echo "step 3"
echo

module purge
module load statistics/R/4.5.0
module load compilers/gcc/15.1.0

Rscript $plot_fst_R $sliding_window_out $diagnostic_region $species $prefix_plot 

