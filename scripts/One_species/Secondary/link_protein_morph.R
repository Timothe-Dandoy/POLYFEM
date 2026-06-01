library("Biostrings")

input_prot_REF <- commandArgs(trailingOnly = TRUE)[1]
input_prot_GUA <- commandArgs(trailingOnly = TRUE)[2]
input_pheno    <- commandArgs(trailingOnly = TRUE)[3]
folder_out     <- commandArgs(trailingOnly = TRUE)[4]

# 1) Link Phenotype - Individuals
Phenotypes <- read.table( 
  input_pheno,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)

colnames(Phenotypes) <- c("ID", "Phenotype")

IDs = Phenotypes$ID

#Phenotypes$Phenotype[Phenotypes$ID == "T14419"] = NA           # Remove of the three individuals with missmatch phenotype in the haplotype network
#Phenotypes$Phenotype[Phenotypes$ID == "T14457"] = NA
#Phenotypes$Phenotype[Phenotypes$ID == "T14432"] = NA

Ligne = Phenotypes[Phenotypes$Phenotype == "S" & !is.na(Phenotypes$Phenotype),]
Male  = Phenotypes[Phenotypes$Phenotype == "RB" & !is.na(Phenotypes$Phenotype),]

All_Haplotypes = list()

# 2) Get the ednrb sequences

seqs <- readAAStringSet(input_prot_GUA)

Sequences <- as.character(seqs)

seq_list  <- strsplit(Sequences, split = "")
names(seq_list) <- sub("_.*", "", names(seq_list))

seq_lengths <- sapply(seq_list, length)
max_length  <- max(seq_lengths)


# 3) Link Phenotype - Haplotype

Maj_L <- c()
Maj_M <- c()

All_L <- c()
All_M <- c()

for (pos in 1:max_length){
    AA_L = c()
    for (el in Ligne$ID){
        seq = seq_list[[el]]
        if (length(seq) < pos  ){
            AA_L = c(AA_L,NA)
        } else{
            AA = seq[pos]
            AA_L = c(AA_L,AA)
        }
    }
    AA_M = c()
    for (el in Male$ID){
        seq = seq_list[[el]]
        if (length(seq) < pos  ){
            AA_M = c(AA_M,NA)
        } else{
            AA = seq[pos]
            AA_M = c(AA_M,AA)
        }
    }
    
    # Get the predominant amino acid at the position, for each morph
    Maj_L[pos] <- names(sort(table(AA_L), decreasing = TRUE))[1]
    Maj_M[pos] <- names(sort(table(AA_M), decreasing = TRUE))[1]
  
    # 1 if all individuals of the same morphs have the same amino acid at that position, 0 otherwise
    All_L[pos] <- ifelse(length(table(AA_L[!is.na(AA_L)])) == 1,1,0)
    All_M[pos] <- ifelse(length(table(AA_M[!is.na(AA_M)])) == 1,1,0)

}

# Outputting the results 
OutPut <- data.frame(
  Position = 1:max_length,
  Maj_L    = Maj_L,
  Maj_M    = Maj_M,
  All_L    = All_L,
  All_M    = All_M
)


OutPut_diff <- OutPut[OutPut$Maj_L != OutPut$Maj_M,]                                   # When the majoritary AA of ligne and male are different

OutPut_megadiff <- OutPut_diff[(OutPut_diff$All_L == 1 & OutPut_diff$All_M == 1),]     # When the majoritary AA of ligne and male are different & all lignes have the same AA ,and same for males

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