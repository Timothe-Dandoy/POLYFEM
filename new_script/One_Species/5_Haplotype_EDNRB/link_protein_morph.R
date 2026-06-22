#-----------------------------------------------------------------
# Packages
#-----------------------------------------------------------------

library(Biostrings)

#-----------------------------------------------------------------
# Data
#-----------------------------------------------------------------

input_prot_REF <- commandArgs(trailingOnly = TRUE)[1]
input_prot <- commandArgs(trailingOnly = TRUE)[2]
input_pheno    <- commandArgs(trailingOnly = TRUE)[3]
folder_out     <- commandArgs(trailingOnly = TRUE)[4]

Phenotypes <- read.table( 
  input_pheno,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)


#-----------------------------------------------------------------
# Scripts
#-----------------------------------------------------------------

# 1) Link Phenotype - Individuals

colnames(Phenotypes) <- c("ID", "Phenotype")
IDs = Phenotypes$ID

Striped     = Phenotypes[Phenotypes$Phenotype == "S" & !is.na(Phenotypes$Phenotype),]
Reticulated = Phenotypes[Phenotypes$Phenotype == "RB" & !is.na(Phenotypes$Phenotype),]

All_Haplotypes = list()

# 2) Get the ednrb protein sequences

seqs <- readAAStringSet(input_prot)
Sequences <- as.character(seqs)

seq_list  <- strsplit(Sequences, split = "")
names(seq_list) <- sub("_.*", "", names(seq_list))

seq_lengths <- sapply(seq_list, length)
max_length  <- max(seq_lengths)


# 3) Link Phenotype - Haplotype

Maj_S <- c()
Maj_RB <- c()

All_S <- c()
All_RB <- c()

# Look at all positions
for (pos in 1:max_length){
    # Get the aminoacid at that position for all striped females
    AA_S = c()
    for (el in Striped$ID){
        seq = seq_list[[el]]
        if (length(seq) < pos  ){
            AA_S = c(AA_S,NA)
        } else{
            AA = seq[pos]
            AA_S = c(AA_S,AA)
        }
    }
    # Get the aminoacid at that position for all RB females
    AA_RB = c()
    for (el in Reticulated$ID){
        seq = seq_list[[el]]
        if (length(seq) < pos  ){
            AA_RB = c(AA_RB,NA)
        } else{
            AA = seq[pos]
            AA_RB = c(AA_RB,AA)
        }
    }
    
    # Get the predominant amino acid at the position, for each morph
    Maj_S[pos] <- names(sort(table(AA_S), decreasing = TRUE))[1]
    Maj_RB[pos] <- names(sort(table(AA_RB), decreasing = TRUE))[1]
  
    # 1 if all individuals of the same morphs have the same amino acid at that position, 0 otherwise
    All_S[pos] <- ifelse(length(table(AA_S[!is.na(AA_S)])) == 1,1,0)
    All_RB[pos] <- ifelse(length(table(AA_RB[!is.na(AA_RB)])) == 1,1,0)

}

# Outputting the results 
OutPut <- data.frame(
  Position = 1:max_length,
  Maj_S    = Maj_S,
  Maj_RB    = Maj_RB,
  All_S    = All_S,
  All_RB    = All_RB
)


OutPut_diff <- OutPut[OutPut$Maj_S != OutPut$Maj_RB,]   # When the majoritary AA of S and RB are different

OutPut_megadiff <- OutPut_diff[(OutPut_diff$All_S == 1 & OutPut_diff$All_RB == 1),]     # When the majoritary AA of S and RB are different & all S have the same AA ,and same for RB

write.table(
  OutPut,
  file = paste0(folder_out,"_OutPut_prot.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  OutPut_diff,
  file = paste0(folder_out,"_OutPut_prot_diff.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  OutPut_megadiff,
  file = paste0(folder_out,"_OutPut_prot_megadiff.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)