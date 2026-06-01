library(readxl)
library(FactoMineR)
library(factoextra)

input_pheno     <- commandArgs(trailingOnly = TRUE)[1]
species         <- commandArgs(trailingOnly = TRUE)[2]
bam_list        <- commandArgs(trailingOnly = TRUE)[3]
bam_list_S      <- commandArgs(trailingOnly = TRUE)[4]
bam_list_RB     <- commandArgs(trailingOnly = TRUE)[5]
individual_file <- commandArgs(trailingOnly = TRUE)[6]
ind_pheno_file  <- commandArgs(trailingOnly = TRUE)[7]
pop_file        <- commandArgs(trailingOnly = TRUE)[8]
pheno_file      <- commandArgs(trailingOnly = TRUE)[9]

# Get data
BD = read_xlsx(path = input_pheno)
BD <- BD[BD$Species == species,]
BD <- BD[!is.na(BD$'genetic data') & BD$'genetic data' == "Y",]

IDs = BD$Sample_ID
Pheno_quali = BD$Phenotype
IDs_S = IDs[!is.na(Pheno_quali) & Pheno_quali  == "S"]
IDs_RB = IDs[!is.na(Pheno_quali) & Pheno_quali  == "RB"]
pops = BD$Population
pops_uni = unique(pops)
pops = match(pops,pops_uni)
bd_pheno <- BD[,c(1,4)]

# PCA
BD$`F - spot`            <- as.numeric(BD$`F - spot`)
BD$`DL - width`          <- as.numeric(BD$`DL - width`)
BD$`DL - discontinuity`  <- as.numeric(BD$`DL - discontinuity`)
BD$`SDS - width`         <- as.numeric(BD$`SDS - width`)
BD$`SDS - discontinuity` <- as.numeric(BD$`SDS - discontinuity`)
BD$`SDS - spot`          <- as.numeric(BD$`SDS - spot`)

BD$`LF - white line` <- ifelse(BD$`LF - white line` == "C", 0,
                                  ifelse(BD$`LF - white line` == "D", 1, NA))

BD$`SDS - regular border` <- ifelse(BD$`SDS - regular border` == "Y", 0,
                                  ifelse(BD$`SDS - regular border` == "N", 1, NA))

BD$`SDS - isolated white scales` <- ifelse(BD$`SDS - isolated white scales` == "N", 0,
                                  ifelse(BD$`SDS - isolated white scales` == "Y", 1, NA))

Var_ACP <- c("LF - white line","SDS - regular border","SDS - isolated white scales","F - spot","DL - width",
             "DL - discontinuity","SDS - width","SDS - discontinuity", "SDS - spot")

BD_clean <- BD[complete.cases(BD[,Var_ACP]), ]
BD_clean <- BD_clean[BD_clean$Sample_ID != "T13969", ]

resPCA <- PCA(BD_clean[,Var_ACP], scale.unit = TRUE, ncp = 5, graph = FALSE)
coord_axe1 <- resPCA$ind$coord[,1]
print(coord_axe1)

dim1 = c()
count = 1
for (el in IDs){
  if (el %in% BD_clean$Sample_ID){
    dim1 = c(dim1,coord_axe1[count])
    count = count + 1 
  } else {dim1 = c(dim1,-999)}
}

# Outputs
Bam = paste0('../../../../bam_analysis/',IDs,'.bam')
writeLines(Bam, bam_list)
Bam_S = paste0('../../../../bam_analysis/',IDs_S,'.bam')
writeLines(Bam_S, bam_list_S)
Bam_RB = paste0('../../../../bam_analysis/',IDs_RB,'.bam')
writeLines(Bam_RB, bam_list_RB)
writeLines(IDs, individual_file)
writeLines(as.character(pops),pop_file)
writeLines(as.character(dim1),pheno_file)

write.table(
  bd_pheno,
  file = ind_pheno_file,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
