#!/bin/bash -l
#SBATCH -J Haplotype_EDNRB_All_Species
#SBATCH -o ../../logs/haplotype_EDNRB_allspe.log
#SBATCH -e ../../logs/haplotype_EDNRB_allspe.log
#SBATCH -t 24:00:00
#SBATCH -p workq
#SBATCH --cpus-per-task 8
#SBATCH --mem=100G

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variables

#Data
Name_Chr=NC_041330.1    # Chromosome Z
Region=30556153-30574830           # EDNRB gene

#Inputs
tmp_folder_pheno=../../../tmp/All_species/Files_pheno/
tmp_folder_genet=../../../tmp/All_species/Files_genet/
file_pheno=../../../data/Pheno.xlsx                                            
fasta_ref_chrz=../../../data/P_muralis/ref_chrZ.fa
fasta_ref_exon=../../../data/P_muralis/podmur1_0_ednrb_exons.fa
Bam_list=${tmp_folder_pheno}bam_allspecies_female_filtered.fofn 

#Temporary
fasta_haplotype=${tmp_folder_genet}haplotypes_ednrb_allspe.fa
fasta_allind_ednrb=${tmp_folder_genet}allspecies_ednrb_V2.fa
fasta_allind_ednrb_aligned=${tmp_folder_genet}allspecies_ednrb_aligned_V2.fa

out_blast_exon=${tmp_folder_genet}Allind_Exons_ednrb.out

region_file_ex1=${tmp_folder_genet}region_exon_1.txt
region_file_ex8=${tmp_folder_genet}region_exon_8.txt
region_ednrb=${tmp_folder_genet}region_ednrb.txt

#Outputs
nexus_allind_ednrb_aligned=../../../results/All_species/Haplotype_EDNRB/allspecies_ednrb_nexus_V2.nex
nexus_allind_no_pheno=../../../results/All_species/Haplotype_EDNRB/allspecies_ednrb_nexus_V2_no_pheno.nex

#Secondary scripts
Make_RegionFiles_exon_R=../Secondary/make_regionfiles_exons.R
Make_RegionFile_ednrb_R=../Secondary/make_regionfile_ednrb.R
Analyse_R=../Secondary/Fasta_to_Nexus.R

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Get EDNRB sequences of every individuals
echo ""
echo "step 1"
echo ""

module load bioinfo/samtools/1.23

> "$fasta_haplotype"

while read line
do
    File=$line
    Individual=$(basename "$File" .bam)
    samtools consensus -r ${Name_Chr}:${Region} -T "$fasta_ref_chrz" -a -f fasta "$File" | \
    sed "s/^>.*/>$Individual/" >> "$fasta_haplotype"
done < "$Bam_list"

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2. Match the sequences to the reference exons in P. muralis
echo ""
echo "step 2"
echo ""

module purge
module load bioinfo/NCBI_Blast+/2.2.28

# Database creation
makeblastdb -in ${fasta_haplotype} -dbtype nucl      

# blast all individual's sequences in all exons  
blastn -db ${fasta_haplotype} -query ${fasta_ref_exon} -out ${out_blast_exon} -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"   

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3. Build the fastas for EDNRB introns and exons of every individuals
echo ""
echo "step 3"
echo ""

module purge
module load bioinfo/Seqtk/1.3
module load statistics/R/4.5.0 

# Make the region file needed for each exon: ID_ind Exon_start Exon_stop  | same for introns
Rscript $Make_RegionFiles_exon_R $out_blast_exon $tmp_folder_genet                

for exon in 1 2 3 4 5 6 7 8
do    
    region_exon=${tmp_folder_genet}region_exon_${exon}.txt
    exon_allind=${tmp_folder_genet}allind_exon_${exon}.fa
    # Output one fasta per exons, with the sequence of each individual
    seqtk subseq ${fasta_haplotype} ${region_exon} | seqtk seq -r - > ${exon_allind}                
done

for intron in 1 2 3 4 5 6 7
do    
    region_intron=${tmp_folder_genet}region_intron_${intron}.txt
    intron_allind=${tmp_folder_genet}allind_intron_${intron}.fa
     # Output one fasta per intron, with the sequence of each individual
    seqtk subseq ${fasta_haplotype} ${region_intron} | seqtk seq -r - > ${intron_allind}           
done

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 4. EDNRB sequence reconstruction 
echo ""
echo "step 4"
echo ""

module purge
module load statistics/R/4.5.0
module load bioinfo/Seqtk/1.3

Rscript $Make_RegionFile_ednrb_R $region_file_ex1 $region_file_ex8 $region_ednrb

# Output one fasta with the edrnb sequence of each individual
seqtk subseq ${fasta_haplotype} ${region_ednrb} > ${fasta_allind_ednrb}                

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 5. EDNRB sequence alignement 
echo ""
echo "step 5"
echo ""

module purge
module load bioinfo/MAFFT/7.505

mafft --auto ${fasta_allind_ednrb} > ${fasta_allind_ednrb_aligned}  

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 6. Build the nexus file with the EDNRB sequence, the species and the phenotype of each individual
echo ""
echo "step 6"
echo ""

module purge
module load statistics/R/4.5.0
module load compilers/gcc/15.1.0

Rscript $Analyse_R $fasta_allind_ednrb_aligned $nexus_allind_ednrb_aligned $nexus_allind_no_pheno $file_pheno