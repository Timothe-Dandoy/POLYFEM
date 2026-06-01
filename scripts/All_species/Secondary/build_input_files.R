library(readxl)

table_all <- commandArgs(trailingOnly = TRUE)[1]
bam_list  <- commandArgs(trailingOnly = TRUE)[2]

BD = read_xlsx(path =table_all)
BD <- BD[!is.na(BD$'genetic data') & BD$'genetic data' == "Y",]

IDs = BD$Sample_ID

Bam = paste0('../../../../bam_analysis/',IDs,'.bam')

writeLines(Bam, bam_list)
