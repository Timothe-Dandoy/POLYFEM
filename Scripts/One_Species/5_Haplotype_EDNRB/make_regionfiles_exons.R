#-----------------------------------------------------------------
# Data
#-----------------------------------------------------------------

output_exons <- commandArgs(trailingOnly = TRUE)[1]
tmp_folder   <- commandArgs(trailingOnly = TRUE)[2]

Allind_exons <- read.table(
  output_exons,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)

#-----------------------------------------------------------------
# Scripts
#-----------------------------------------------------------------

colnames(Allind_exons) <- c("Exon", "IDs",3:6,"end","start",9,10)

Allind_exons$Exon <- sub(".*Exon_", "", Allind_exons$Exon)

Allind_exons <- Allind_exons[c(1:6,8,7,9,10)]         # Switching the start and the end 

for (ex in 1:8){
    region_exon = paste0(tmp_folder,"_region_exon_",ex,".txt")
    
    Allind_1exon = Allind_exons[Allind_exons$Exon == ex,]
    Allind_1exon = Allind_1exon[c(2,7,8)]
    Allind_1exon[,2] = Allind_1exon[,2]-1

    Allind_1exon <- Allind_1exon[order(Allind_1exon$IDs), ]

    write.table(
        Allind_1exon,
        file = region_exon,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )
}

for (int in 1:7){
    region_intron = paste0(tmp_folder,"_region_intron_",int,".txt")
    
    Allind_1exon_1 = Allind_exons[Allind_exons$Exon == int,]      #looks at the previous exon
    Allind_1exon_2 = Allind_exons[Allind_exons$Exon == int+1,]    #looks at the next exon

    Allind_1exon_1 <- Allind_1exon_1[order(Allind_1exon_1$IDs), ]
    Allind_1exon_2 <- Allind_1exon_2[order(Allind_1exon_2$IDs), ]

    Allind_1intron <- data.frame(
        Ids = Allind_1exon_1$IDs,
        Start = Allind_1exon_2[,8],             # The intron start at the end of the previous exon
        End = Allind_1exon_1[,7] - 1                # and end at the start of the next one.
    )

    write.table(
        Allind_1intron,
        file = region_intron,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )
}


