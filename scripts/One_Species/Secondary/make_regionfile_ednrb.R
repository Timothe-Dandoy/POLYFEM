#-----------------------------------------------------------------
# Data
#-----------------------------------------------------------------

exon_1_file   <- commandArgs(trailingOnly = TRUE)[1]
exon_8_file   <- commandArgs(trailingOnly = TRUE)[2]
region_ednrb  <- commandArgs(trailingOnly = TRUE)[3]

exon_1 <- read.table(
  exon_1_file,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)

exon_8 <- read.table(
  exon_8_file,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)

#-----------------------------------------------------------------
# Scripts
#-----------------------------------------------------------------

colnames(exon_1) <- c("IDs", "start","end")
colnames(exon_8) <- c("IDs", "start","end")

exon_1 <- exon_1[order(exon_1$IDs), ]
exon_8 <- exon_8[order(exon_8$IDs), ]

ednrb = data.frame(
    IDs   = exon_1$IDs,
    start = exon_8$start,
    end   = exon_1$end
)

write.table(
    ednrb,
    file = region_ednrb,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
)


