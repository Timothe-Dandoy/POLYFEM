#-----------------------------------------------------------------
# Packages
#-----------------------------------------------------------------

library(ape)
library(readxl)

#-----------------------------------------------------------------
# Data
#-----------------------------------------------------------------

input_fasta    <- commandArgs(trailingOnly = TRUE)[1]
output_nexus   <- commandArgs(trailingOnly = TRUE)[2]
output_nopheno <- commandArgs(trailingOnly = TRUE)[3]
table_all      <- commandArgs(trailingOnly = TRUE)[4]
 
ednrb <- read.dna(input_fasta, format = "fasta")
rownames(ednrb) <- sub(":.*", "", rownames(ednrb))

Data_base <- read_xlsx(path = table_all)

#-----------------------------------------------------------------
# Scripts
#-----------------------------------------------------------------

# Nexus files creation
# with phenotype
write.nexus.data(ednrb, file = output_nexus, format = "dna",interleaved = FALSE)
#without phenotype
write.nexus.data(ednrb, file = output_nopheno, format = "dna",interleaved = FALSE)

# Collect the phenotype and the species of all individuals
Data_base <- Data_base[!is.na(Data_base$'genetic data') & Data_base$'genetic data' == "Y",]
Data_base <- Data_base[, c("Sample_ID", "Species", "Phenotype")]

phenotype_levels <- c("S","RB")
species_levels   <- c("P_bocagei", "P_carbonelli", "P_guadarramae","P_hispanicus", "P_liolepis", "P_lusitanicus","P_tunesiacus","P_vaucheri", "P_virescens")

Data_base$Ligne = ifelse(is.na(Data_base$Phenotype), "?" , ifelse(Data_base$Phenotype == "S"  , 1 , ifelse(Data_base$Phenotype %in% phenotype_levels, 0 , "?")))
Data_base$Morph = ifelse(is.na(Data_base$Phenotype), "?" , ifelse(Data_base$Phenotype == "RB" , 1 , ifelse(Data_base$Phenotype %in% phenotype_levels, 0 , "?")))

Data_base$bocagei     = ifelse(is.na(Data_base$Species), "?" , ifelse(Data_base$Species == "P_bocagei"    , 1 , ifelse(Data_base$Species %in% species_levels, 0 , "?")))
Data_base$carbonelli  = ifelse(is.na(Data_base$Species), "?" , ifelse(Data_base$Species == "P_carbonelli" , 1 , ifelse(Data_base$Species %in% species_levels, 0 , "?")))
Data_base$guadarramae = ifelse(is.na(Data_base$Species), "?" , ifelse(Data_base$Species == "P_guadarramae", 1 , ifelse(Data_base$Species %in% species_levels, 0 , "?")))
Data_base$hispanicus  = ifelse(is.na(Data_base$Species), "?" , ifelse(Data_base$Species == "P_hispanicus" , 1 , ifelse(Data_base$Species %in% species_levels, 0 , "?")))
Data_base$liolepis    = ifelse(is.na(Data_base$Species), "?" , ifelse(Data_base$Species == "P_liolepis"   , 1 , ifelse(Data_base$Species %in% species_levels, 0 , "?")))
Data_base$lusitanicus = ifelse(is.na(Data_base$Species), "?" , ifelse(Data_base$Species == "P_lusitanicus", 1 , ifelse(Data_base$Species %in% species_levels, 0 , "?")))
Data_base$tunesiacus  = ifelse(is.na(Data_base$Species), "?" , ifelse(Data_base$Species == "P_tunesiacus" , 1 , ifelse(Data_base$Species %in% species_levels, 0 , "?")))
Data_base$vaucheri    = ifelse(is.na(Data_base$Species), "?" , ifelse(Data_base$Species == "P_vaucheri"   , 1 , ifelse(Data_base$Species %in% species_levels, 0 , "?")))
Data_base$virescens   = ifelse(is.na(Data_base$Species), "?" , ifelse(Data_base$Species == "P_virescens"  , 1 , ifelse(Data_base$Species %in% species_levels, 0 , "?")))

# Adding the covariates to the nexus files
Traits_matrix = cbind(Data_base$Sample_ID,
                      Data_base$Ligne,Data_base$Morph,
                      Data_base$bocagei,Data_base$carbonelli,Data_base$guadarramae,Data_base$hispanicus,Data_base$liolepis,Data_base$lusitanicus,Data_base$tunesiacus,Data_base$vaucheri,Data_base$virescens
                      )

Text_matrix = apply(Traits_matrix, 1, function(x) paste(x, collapse=" "))

traits_block <- c(
  "BEGIN TRAITS;",
  "Dimensions NTRAITS=11;",
  "Format labels=yes missing=? separator=spaces;",
  "TraitLabels S RB Bocagei Carbonelli Guadarramae Hispanicus Liolepis Lusitanicus Tunesiacus Vaucheri Virescens;",
  "Matrix",
  paste0("\t",Text_matrix),
  ";",
  "END;"
)

write(traits_block,
      file = output_nexus,
      append = TRUE) 

Traits_matrix_2 = cbind(Data_base$Sample_ID,
                      Data_base$bocagei,Data_base$carbonelli,Data_base$guadarramae,Data_base$hispanicus,Data_base$liolepis,Data_base$lusitanicus,Data_base$tunesiacus,Data_base$vaucheri,Data_base$virescens
                      )

Text_matrix_2= apply(Traits_matrix_2, 1, function(x) paste(x, collapse=" "))

traits_block_2 <- c(
  "BEGIN TRAITS;",
  "Dimensions NTRAITS=9;",
  "Format labels=yes missing=? separator=spaces;",
  "TraitLabels Bocagei Carbonelli Guadarramae Hispanicus Liolepis Lusitanicus Tunesiacus Vaucheri Virescens;",
  "Matrix",
  paste0("\t",Text_matrix_2),
  ";",
  "END;"
)

write(traits_block_2,
      file = output_nopheno,
      append = TRUE)
