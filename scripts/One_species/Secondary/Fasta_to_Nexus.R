library(ape)
library(readxl)

input_fasta         <- commandArgs(trailingOnly = TRUE)[1]
output_nexus        <- commandArgs(trailingOnly = TRUE)[2]
Pheno_ind           <- commandArgs(trailingOnly = TRUE)[3]

# Nexus file creation 
ednrb <- read.dna(input_fasta, format = "fasta")
rownames(ednrb) <- sub(":.*", "", rownames(ednrb))
write.nexus.data(ednrb, file = output_nexus, format = "dna",interleaved = FALSE)

# Adding the covariate => phenotype  

Phenotypes <- read.table(
  Pheno_ind,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)

colnames(Phenotypes) <- c("ID", "Phenotype")

phenotype_levels <- c("S","RB")

Phenotypes$Ligne = ifelse(is.na(Phenotypes$Phenotype), "?" , ifelse(Phenotypes$Phenotype == "S" , 1 , ifelse(Phenotypes$Phenotype %in% phenotype_levels, 0 , "?")))
Phenotypes$Male  = ifelse(is.na(Phenotypes$Phenotype), "?" , ifelse(Phenotypes$Phenotype == "RB" , 1 , ifelse(Phenotypes$Phenotype %in% phenotype_levels, 0 , "?")))

Traits_matrix = cbind(Phenotypes$ID,
                      Phenotypes$Ligne,Phenotypes$Male
                      )

Text_matrix = apply(Traits_matrix, 1, function(x) paste(x, collapse=" "))

traits_block <- c(
  "BEGIN TRAITS;",
  "Dimensions NTRAITS=2;",
  "Format labels=yes missing=? separator=spaces;",
  "TraitLabels Lignee Male;",
  "Matrix",
  paste0("\t",Text_matrix),
  ";",
  "END;"
)

write(traits_block,
      file = output_nexus,
      append = TRUE)


   