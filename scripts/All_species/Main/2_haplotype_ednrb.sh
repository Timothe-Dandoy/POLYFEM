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

Name_Chr=NC_041330.1
Region=30556153-30574830
tmp_folder_pheno=../../../tmp/All_species/Files_pheno/
tmp_folder_genet=../../../tmp/All_species/Files_genet/

file_pheno=../../../data/Pheno.xlsx                                            
fasta_ref_chrz=../../../data/P_muralis/ref_chrZ.fa
fasta_ref_exon=../../../data/P_muralis/podmur1_0_ednrb_exons.fa
Bam_list=${tmp_folder_pheno}bam_allspecies_female_filtered.fofn 

fasta_haplotype=${tmp_folder_genet}haplotypes_ednrb_allspe.fa
out_blast_exon=${tmp_folder_genet}Allind_Exons_ednrb.out
fasta_allind_ednrb=${tmp_folder_genet}allspecies_ednrb_V2.fa
fasta_allind_ednrb_aligned=${tmp_folder_genet}allspecies_ednrb_aligned_V2.fa
nexus_allind_ednrb_aligned=../../../results/All_species/Haplotype_EDNRB/allspecies_ednrb_nexus_V2.nex
nexus_allind_no_pheno=../../../results/All_species/Haplotype_EDNRB/allspecies_ednrb_nexus_V2_no_pheno.nex

region_file_ex1=${tmp_folder_genet}region_exon_1.txt
region_file_ex8=${tmp_folder_genet}region_exon_8.txt
region_ednrb=${tmp_folder_genet}region_ednrb.txt

Make_RegionFiles_exon_R=../Secondary/make_regionfiles_exons.R
Make_RegionFile_ednrb_R=../Secondary/make_regionfile_ednrb.R
Analyse_R=../Secondary/Fasta_to_Nexus.R

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Get Haplotype
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
# 3. Match Haplotypes to Reference exons
echo ""
echo "step 3"
echo ""

module purge
module load bioinfo/NCBI_Blast+/2.2.28

makeblastdb -in ${fasta_haplotype} -dbtype nucl      # Database creation
blastn -db ${fasta_haplotype} -query ${fasta_ref_exon} -out ${out_blast_exon} -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"   # blast all individual's sequences in all exons  

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 4. Fasta for Introns and Exons
echo ""
echo "step 4"
echo ""

module purge
module load bioinfo/Seqtk/1.3
module load statistics/R/4.5.0 

Rscript $Make_RegionFiles_exon_R $out_blast_exon $tmp_folder_genet                # Make the region file needed for each exon: ID_ind Exon_start Exon_stop  | same for introns

for exon in 1 2 3 4 5 6 7 8
do    
    region_exon=${tmp_folder_genet}region_exon_${exon}.txt
    exon_allind=${tmp_folder_genet}allind_exon_${exon}.fa
    seqtk subseq ${fasta_haplotype} ${region_exon} | seqtk seq -r - > ${exon_allind}                # Output one fasta per exons, with the sequence of each individual
done

for intron in 1 2 3 4 5 6 7
do    
    region_intron=${tmp_folder_genet}region_intron_${intron}.txt
    intron_allind=${tmp_folder_genet}allind_intron_${intron}.fa
    seqtk subseq ${fasta_haplotype} ${region_intron} | seqtk seq -r - > ${intron_allind}            # Output one fasta per intron, with the sequence of each individual
done

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 5. EDNRB sequence reconstruction 
echo ""
echo "step 5"
echo ""

module purge
module load statistics/R/4.5.0
module load bioinfo/Seqtk/1.3

Rscript $Make_RegionFile_ednrb_R $region_file_ex1 $region_file_ex8 $region_ednrb
 
seqtk subseq ${fasta_haplotype} ${region_ednrb} > ${fasta_allind_ednrb}                # Output one fasta with the edrnb sequence of each individual

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 6. EDNRB sequence alignement 
echo ""
echo "step 6"
echo ""

module purge
module load bioinfo/MAFFT/7.505

mafft --auto ${fasta_allind_ednrb} > ${fasta_allind_ednrb_aligned}  

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 7. Fasta to Nexus
echo ""
echo "step 7"
echo ""

module purge
module load statistics/R/4.5.0
module load compilers/gcc/15.1.0

Rscript $Analyse_R $fasta_allind_ednrb_aligned $nexus_allind_ednrb_aligned $nexus_allind_no_pheno $file_pheno