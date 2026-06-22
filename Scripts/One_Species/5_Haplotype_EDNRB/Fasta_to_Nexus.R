#-----------------------------------------------------------------
# Packages
#-----------------------------------------------------------------

library(ape)
library(readxl)

#-----------------------------------------------------------------
# Data
#-----------------------------------------------------------------

input_fasta         <- commandArgs(trailingOnly = TRUE)[1]
output_nexus        <- commandArgs(trailingOnly = TRUE)[2]
Pheno_ind           <- commandArgs(trailingOnly = TRUE)[3]

ednrb <- read.dna(input_fasta, format = "fasta")
rownames(ednrb) <- sub(":.*", "", rownames(ednrb))

Phenotypes <- read.table(
  Pheno_ind,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)

#-----------------------------------------------------------------
# Scripts
#-----------------------------------------------------------------

# Nexus file creation
write.nexus.data(ednrb, file = output_nexus, format = "dna",interleaved = FALSE)

# Collect the phenotype of all individuals
colnames(Phenotypes) <- c("ID", "Phenotype")

phenotype_levels <- c("S","RB")

Phenotypes$S   = ifelse(is.na(Phenotypes$Phenotype), "?" , ifelse(Phenotypes$Phenotype == "S" , 1  , ifelse(Phenotypes$Phenotype %in% phenotype_levels, 0 , "?")))
Phenotypes$RB  = ifelse(is.na(Phenotypes$Phenotype), "?" , ifelse(Phenotypes$Phenotype == "RB" , 1 , ifelse(Phenotypes$Phenotype %in% phenotype_levels, 0 , "?")))

# Adding the covariate to the nexus file
Traits_matrix = cbind(Phenotypes$ID,
                      Phenotypes$S,
                      Phenotypes$RB)

Text_matrix = apply(Traits_matrix, 1, function(x) paste(x, collapse=" "))

traits_block <- c(
  "BEGIN TRAITS;",
  "Dimensions NTRAITS=2;",
  "Format labels=yes missing=? separator=spaces;",
  "TraitLabels Striped ReticulatedBlotched;",
  "Matrix",
  paste0("\t",Text_matrix),
  ";",
  "END;"
)

write(traits_block,
      file = output_nexus,
      append = TRUE)


   