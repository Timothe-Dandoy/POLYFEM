#!/bin/bash -l
#SBATCH -J Haplotype_EDNRB
#SBATCH -o ../../logs/haplotype_EDNRB.log
#SBATCH -e ../../logs/haplotype_EDNRB.log
#SBATCH -t 24:00:00
#SBATCH -p workq
#SBATCH --cpus-per-task 8
#SBATCH --mem=100G

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Variables

species=$1
nickname=$(echo ${species} | cut -d'_' -f2 | cut -c1-3) 
Name_Chr=NC_041330.1
Region=30556153-30574830

data_folder=../../../data/${species}/
tmp_folder_pheno=../../../tmp/${species}/Files_pheno/${nickname}
tmp_folder_genet=../../../tmp/${species}/Files_genet/${nickname}
results_folder=../../../results/${species}/Haplotype_EDNRB/${nickname}

Bam_list=${tmp_folder_pheno}_bam_list_filtered.fofn
individual_file=${tmp_folder_pheno}_individual.txt
pheno_ind=${tmp_folder_pheno}_individual_phenotype.txt

fasta_ref_chrz=../../../data/P_muralis/ref_chrZ.fa
fasta_ref_exon=../../../data/P_muralis/podmur1_0_ednrb_exons.fa
fasta_ref_coding=../../../data/P_muralis/podmur1_0_ednrb_coding.fa
fasta_ref_protein=../../../data/P_muralis/Ref_protein_ednrb_PodMur1_0.fa

Prefix_exon=${tmp_folder_genet}_allind_exon_
Prefix_intron=${tmp_folder_genet}_allind_intron_
Suffix=_aligned.fa

fasta_haplotype=${tmp_folder_genet}_haplotypes_allind.fa
fasta_allind_mRNA=${tmp_folder_genet}_allind_mRNA.fa
fasta_allind_coding=${tmp_folder_genet}_allind_coding_ednrb.fa
fasta_allind_protein=${tmp_folder_genet}_allind_protein_ednrb.fa
fasta_allind_ednrb=${tmp_folder_genet}_allind_ednrb.fa
fasta_allind_ednrb_aligned=${tmp_folder_genet}_allind_ednrb_aligned.fa
nexus_allind_ednrb_aligned=${results_folder}_ednrb_nexus.nex

out_blast_exon=${tmp_folder_genet}_allind_exons_ednrb.out
out_blast_coding=${tmp_folder_genet}_allind_coding_ednrb.out

region_coding=${tmp_folder_genet}_region_coding.txt
region_file_ex1=${tmp_folder_genet}_region_exon_1.txt
region_file_ex8=${tmp_folder_genet}_region_exon_8.txt
region_ednrb=${tmp_folder_genet}_region_ednrb.txt

Make_RegionFiles_exon_R=../Secondary/make_regionfiles_exons.R
Make_RegionFile_coding_R=../Secondary/make_regionfile_coding.R
Make_RegionFile_ednrb_R=../Secondary/make_regionfile_ednrb.R
prot_to_morph_R=../Secondary/link_protein_morph.R
nuc_to_morph_R=../Secondary/link_nucleotide_morph.R
fasta_to_nexus_R=../Secondary/Fasta_to_Nexus.R

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Get Haplotypes
echo "step 1"

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
# 2. Match Haplotypes to Reference exons
echo "step 2"

module purge
module load bioinfo/NCBI_Blast+/2.2.28

makeblastdb -in ${fasta_haplotype} -dbtype nucl      # Database creation
blastn -db ${fasta_haplotype} -query ${fasta_ref_exon} -out ${out_blast_exon} -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"   # blast all individual's sequences in all exons  

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3. Fasta for Introns and Exons
echo "step 3"

module purge
module load bioinfo/Seqtk/1.3
module load statistics/R/4.5.0
module load compilers/gcc/15.1.0

Rscript $Make_RegionFiles_exon_R $out_blast_exon $tmp_folder_genet                # Make the region file needed for each exon: ID_ind Exon_start Exon_stop  | same for introns

for exon in 1 2 3 4 5 6 7 8 
do    
    region_exon=${tmp_folder_genet}_region_exon_${exon}.txt
    exon_allind=${tmp_folder_genet}_allind_exon_${exon}.fa
    seqtk subseq ${fasta_haplotype} ${region_exon} | seqtk seq -r - > ${exon_allind}                # Output one fasta per exons, with the sequence of each individual
done

for intron in 1 2 3 4 5 6 7
do    
    region_intron=${tmp_folder_genet}_region_intron_${intron}.txt
    intron_allind=${tmp_folder_genet}_allind_intron_${intron}.fa
    seqtk subseq ${fasta_haplotype} ${region_intron} | seqtk seq -r - > ${intron_allind}            # Output one fasta per intron, with the sequence of each individual
done

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 4. One fasta for all Exons -> mRNA
echo "step 4"

awk '{print $1}' ${region_file_ex1} | sort -u > ${individual_file}        # Make a file with the list of all individuals IDs

> ${fasta_allind_mRNA}

while read individual
do
    echo ">${individual}" >> ${fasta_allind_mRNA}                      
    
    concatenated_seq=""
    for exon in 1 2 3 4 5 6 7 8
    do
        exon_file=${tmp_folder_genet}_allind_exon_${exon}.fa
        seq=$(grep -A1 "^>${individual}" ${exon_file} | tail -n1)
        concatenated_seq="${concatenated_seq}${seq}"
    done
    echo "${concatenated_seq}" >> ${fasta_allind_mRNA}
done < ${individual_file} 

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 5. Fasta for the coding sequence  
echo "step 5"

module purge
module load bioinfo/NCBI_Blast+/2.2.28
module load bioinfo/Seqtk/1.3
module load statistics/R/4.5.0
module load compilers/gcc/15.1.0

makeblastdb -in ${fasta_allind_mRNA} -dbtype nucl      # Database creation
blastn -db ${fasta_allind_mRNA} -query ${fasta_ref_coding} -out ${out_blast_coding} -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"   # blast all individual's sequences in the coding sequence  

Rscript $Make_RegionFile_coding_R $out_blast_coding $region_coding
 
seqtk subseq ${fasta_allind_mRNA} ${region_coding} > ${fasta_allind_coding}                # Output one fasta with the coding sequence of each individual

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 6. Protein Translation 
echo "step 6"

module purge
module load bioinfo/EMBOSS/6.6.0

# For GUA individuals
sed -i 's/:/_/g' ${fasta_allind_coding}     # change ":" to "_"

transeq -sequence ${fasta_allind_coding} \
        -outseq ${fasta_allind_protein} \
        -frame 1 \
        -table 1

# For the reference P. muralis PodMur 1.0
sed -i 's/:/_/g' ${fasta_ref_coding}     # change ":" to "_"

transeq -sequence ${fasta_ref_coding} \
        -outseq ${fasta_ref_protein} \
        -frame 1 \
        -table 1

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 7. Protein comparison between phenotypes 
echo "step 7"

module purge
module load statistics/R/4.5.0
module load compilers/gcc/15.1.0

Rscript $prot_to_morph_R $fasta_ref_protein $fasta_allind_protein $pheno_ind $results_folder

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 8. Sequences alignment 
echo "step 8"

module purge
module load bioinfo/MAFFT/7.505

for exon in 1 2 3 4 5 6 7 8
do    
    fasta=${tmp_folder_genet}_allind_exon_${exon}.fa
    fasta_aligned=${tmp_folder_genet}_allind_exon_${exon}_aligned.fa
    mafft --auto ${fasta} > ${fasta_aligned}             
done

for intron in 1 2 3 4 5 6 7
do    
    fasta=${tmp_folder_genet}_allind_intron_${intron}.fa
    fasta_aligned=${tmp_folder_genet}_allind_intron_${intron}_aligned.fa
    mafft --auto ${fasta} > ${fasta_aligned}       
done

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 9. Nucleotide comparison between phenotypes 
echo "step 9"

module purge
module load statistics/R/4.5.0
module load compilers/gcc/15.1.0

Rscript $nuc_to_morph_R $Prefix_exon $Prefix_intron $Suffix $pheno_ind $results_folder

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 10. EDNRB sequence reconstruction 
echo "step 10"

module purge
module load statistics/R/4.5.0
module load compilers/gcc/15.1.0
module load bioinfo/Seqtk/1.3

Rscript $Make_RegionFile_ednrb_R $region_file_ex1 $region_file_ex8 $region_ednrb
 
seqtk subseq ${fasta_haplotype} ${region_ednrb} > ${fasta_allind_ednrb}                # Output one fasta with the edrnb sequence of each individual

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 11. EDNRB sequence alignement 
echo "step 11"

module purge
module load bioinfo/MAFFT/7.505

mafft --auto ${fasta_allind_ednrb} > ${fasta_allind_ednrb_aligned}  

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 12. Fasta to Nexus
echo "step 12"

module purge
module load statistics/R/4.5.0
module load compilers/gcc/15.1.0

Rscript $fasta_to_nexus_R $fasta_allind_ednrb_aligned $nexus_allind_ednrb_aligned $pheno_ind