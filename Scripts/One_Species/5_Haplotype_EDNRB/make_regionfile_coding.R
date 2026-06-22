#-----------------------------------------------------------------
# Data
#-----------------------------------------------------------------

output_coding <- commandArgs(trailingOnly = TRUE)[1]
region_coding <- commandArgs(trailingOnly = TRUE)[2]

Allind_coding <- read.table(
  output_coding,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)

#-----------------------------------------------------------------
# Scripts
#-----------------------------------------------------------------

colnames(Allind_coding) <- c("Coding", "IDs",3:10)
    
Allind_coding = Allind_coding[c(2,7,8)]
Allind_coding[,2] = Allind_coding[,2]-1

Allind_coding <- Allind_coding[order(Allind_coding$IDs), ]

write.table(
    Allind_coding,
    file = region_coding,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
)


