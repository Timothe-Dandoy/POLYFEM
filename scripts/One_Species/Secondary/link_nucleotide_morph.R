#-----------------------------------------------------------------
# Packages
#-----------------------------------------------------------------

library(ggplot2)
library(Biostrings)
library(patchwork)

#-----------------------------------------------------------------
# Data
#-----------------------------------------------------------------

Prefix_exon    <- commandArgs(trailingOnly = TRUE)[1]
Prefix_intron  <- commandArgs(trailingOnly = TRUE)[2]
Suffix         <- commandArgs(trailingOnly = TRUE)[3]
in_pheno       <- commandArgs(trailingOnly = TRUE)[4]
folder_out     <- commandArgs(trailingOnly = TRUE)[5]

Phenotypes <- read.table(
  in_pheno,
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

# Remove the three individuals with missmatch phenotype in the haplotype network
Phenotypes$Phenotype[Phenotypes$ID == "T14419"] = NA          
Phenotypes$Phenotype[Phenotypes$ID == "T14457"] = NA
Phenotypes$Phenotype[Phenotypes$ID == "T14432"] = NA

Striped      = Phenotypes[Phenotypes$Phenotype == "S"  & !is.na(Phenotypes$Phenotype),]      
Reticulated  = Phenotypes[Phenotypes$Phenotype == "RB" & !is.na(Phenotypes$Phenotype),]

# 2) Link Phenotype - Haplotype for all exons and introns 
List_out          = list()           # All the positions
List_out_diff     = list()           # Positions with different alleles for S and RB individuals, in majority 
List_out_megadiff = list()           # Positions with different alleles for S and RB individuals, in total


# look intron by intron, exon by exon : i even => intron  |  i odd => exon 
for (i in 1:15){
    # Collect the nucleotide sequence
    Fasta = ifelse(i %% 2 == 1,paste0(Prefix_exon,i %/% 2 + 1,Suffix),paste0(Prefix_intron,i %/% 2,Suffix))
    seqs <- readDNAStringSet(Fasta)
    Sequences <- as.character(seqs)
    seq_list  <- strsplit(Sequences, split = "")
    names(seq_list) <- sub(":.*", "", names(seq_list))
   
    seq_lengths <- sapply(seq_list, length)
    max_length  <- max(seq_lengths)
    
    Maj_S  <- c() 
    Maj_RB <- c()

    All_S  <- c()
    All_RB <- c()

    Nb_A_S    <- c()
    Nb_C_S    <- c()
    Nb_G_S    <- c()
    Nb_T_S    <- c()
    Nb_N_S    <- c()
    Nb_del_S  <- c()

    Nb_A_RB   <- c()
    Nb_C_RB   <- c()
    Nb_G_RB   <- c()
    Nb_T_RB   <- c()
    Nb_N_RB   <- c()
    Nb_del_RB <- c()

    # Look at all positions
    for (pos in 1:max_length){
        # Get the nucleotides at that position for all striped females
        Nucs_S = c()
        for (el in Striped$ID){
            seq = seq_list[[el]]
            if (length(seq) < pos){
                Nucs_S = c(Nucs_S,NA)
            } else{
                Nuc = seq[pos]
                Nucs_S = c(Nucs_S,Nuc)
            }
        }
        # Get the nucleotides at that position for all reticulated blotched females
        Nucs_RB = c()
        for (el in Reticulated$ID){
            seq = seq_list[[el]]
            if (length(seq) < pos  ){
                Nucs_RB = c(Nucs_RB,NA)
            } else{
                Nuc = seq[pos]
                Nucs_RB = c(Nucs_RB,Nuc)
            }
        }
        
        # Get the predominant nucleotide at the position, for each phenotype
        Maj_S[pos]     <- names(sort(table(Nucs_S), decreasing = TRUE))[1]
        Maj_RB[pos]    <- names(sort(table(Nucs_RB), decreasing = TRUE))[1]
    
        # 1 if all individuals of the same phenotype have the same nucleotide at that position, 0 otherwise
        All_S[pos]     <- ifelse(length(table(Nucs_S[!is.na(Nucs_S)])) == 1,1,0)
        All_RB[pos]    <- ifelse(length(table(Nucs_RB[!is.na(Nucs_RB)])) == 1,1,0)

        # Count the number of each nucleotide at the position, for each phenotype
        Nb_A_S[pos]    <- ifelse(is.na(as.numeric(table(Nucs_S)["A"])),0,as.numeric(table(Nucs_S)["A"]))
        Nb_C_S[pos]    <- ifelse(is.na(as.numeric(table(Nucs_S)["C"])),0,as.numeric(table(Nucs_S)["C"]))
        Nb_G_S[pos]    <- ifelse(is.na(as.numeric(table(Nucs_S)["G"])),0,as.numeric(table(Nucs_S)["G"]))
        Nb_T_S[pos]    <- ifelse(is.na(as.numeric(table(Nucs_S)["T"])),0,as.numeric(table(Nucs_S)["T"]))
        Nb_N_S[pos]    <- ifelse(is.na(as.numeric(table(Nucs_S)["N"])),0,as.numeric(table(Nucs_S)["N"]))
        Nb_del_S[pos]  <- ifelse(is.na(as.numeric(table(Nucs_S)["-"])),0,as.numeric(table(Nucs_S)["-"]))

        Nb_A_RB[pos]   <- ifelse(is.na(as.numeric(table(Nucs_RB)["A"])),0,as.numeric(table(Nucs_RB)["A"]))
        Nb_C_RB[pos]   <- ifelse(is.na(as.numeric(table(Nucs_RB)["C"])),0,as.numeric(table(Nucs_RB)["C"]))
        Nb_G_RB[pos]   <- ifelse(is.na(as.numeric(table(Nucs_RB)["G"])),0,as.numeric(table(Nucs_RB)["G"]))
        Nb_T_RB[pos]   <- ifelse(is.na(as.numeric(table(Nucs_RB)["T"])),0,as.numeric(table(Nucs_RB)["T"]))
        Nb_N_RB[pos]   <- ifelse(is.na(as.numeric(table(Nucs_RB)["N"])),0,as.numeric(table(Nucs_RB)["N"]))
        Nb_del_RB[pos] <- ifelse(is.na(as.numeric(table(Nucs_RB)["-"])),0,as.numeric(table(Nucs_RB)["-"]))
    }
    
    out = ifelse(i %% 2 == 1,paste0("ex",i %/% 2 + 1),paste0("int",i %/% 2))

    OutPut <- data.frame(
    Position = 1:max_length,
    Ex_Int   = rep(out, max_length),
    Maj_S    = Maj_S,
    Maj_RB    = Maj_RB,
    All_S    = All_S,
    All_RB    = All_RB,
    Nb_A     = paste(Nb_A_S,"-",Nb_A_RB),
    Nb_C     = paste(Nb_C_S,"-",Nb_C_RB),
    Nb_G     = paste(Nb_G_S,"-",Nb_G_RB),
    Nb_T     = paste(Nb_T_S,"-",Nb_T_RB),
    Nb_N     = paste(Nb_N_S,"-",Nb_N_RB),
    Nb_del   = paste(Nb_del_S,"-",Nb_del_RB)
    )

    List_out[[i]] <- OutPut
}


# 3) Regroup all exons and introns in one sequence

last_pos = 0

border = list()

for (i in 15:1){
    List_out[[i]]$Position = List_out[[i]]$Position + last_pos
    last_pos  = max(List_out[[i]]$Position)
    first_pos = min(List_out[[i]]$Position)
    border[[i]] = c(first_pos,last_pos)
}


# 4) Output the results
OutPut          <- do.call(rbind, List_out)
OutPut          <- OutPut[order(OutPut$Position), ]
OutPut_diff     <- OutPut[OutPut$Maj_S != OutPut$Maj_RB,]                               # When the majoritary alleles of S and RB are different
OutPut_megadiff <- OutPut_diff[(OutPut_diff$All_S == 1 & OutPut_diff$All_RB == 1),]     # When the majoritary alleles of S and RB are different & all S have the same nucleotide ,and same for RB


write.table(
OutPut,
file = paste0(folder_out,"_ednrb_haplotype.txt"),
sep = "\t",
row.names = FALSE,
quote = FALSE
)

write.table(
OutPut_diff,
file = paste0(folder_out,"_ednrb_haplotype_diff.txt"),
sep = "\t",
row.names = FALSE,
quote = FALSE
)

write.table(
OutPut_megadiff,
file = paste0(folder_out,"_ednrb_haplotype_megadiff.txt"),
sep = "\t",
row.names = FALSE,
quote = FALSE
)

# 5) EDNRB gene representation

start_gene = 30556153

EDNRB <- data.frame(
    start = sapply(border, min) + start_gene -1,
    end   = sapply(border, max) + start_gene -1,
    type  = ifelse((1:15)%% 2 == 1,"Exon","Intron"),
    num   = ((1:15)+1) %/% 2
    )

EDNRB$center  <- (EDNRB$start + EDNRB$end) / 2
Pos_Sign      <- data.frame(pos = OutPut_megadiff$Position + start_gene -1)
Region_Sign   <- c(30562000,30556900)

# Plot of the whole EDNRB sequence
p1 <- ggplot(EDNRB) +
  
  geom_segment(aes(x = start, xend = end,
                   y = 0, yend = 0,
                   color = type),
               size = 6) +
  
  geom_segment(data = Pos_Sign,
               aes(x = pos, xend = pos,
                   y = -0.009, yend = 0.009,
                   color = "Diagnostic positions"),
               linewidth = 0.8) +
  
  geom_rect(aes(xmin = Region_Sign[1],
                xmax = Region_Sign[2],
                ymin = -0.022,
                ymax = 0.02,
                color = "Significant region"),
            fill = NA,
            linewidth = 1,
            inherit.aes = FALSE) +

  scale_x_reverse(breaks = c(30557000, 30560000, 30563000, 30566000, 30569000, 30572000, 30575000)) +

  scale_color_manual(values = c(
    "Exon" = "green3",
    "Intron" = "grey",
    "Diagnostic positions" = "red",
    "Significant region" = "blue"
  )) +

  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14)) +
  
  labs(x = "Position",
       color = "Type") +

  coord_cartesian(ylim = c(-0.2, 0.2))

ggsave(
       paste0(folder_out,"_Gene_EDNRB.png"),
       plot = p1 ,
       width = 10,
       height = 5,
       dpi = 150   
   )

# Zoom in the diagnostic region
p2 <- ggplot(EDNRB) +
  
  geom_segment(aes(x = start, xend = end,
                   y = 0, yend = 0,
                   color = type),
               size = 6) +

  geom_segment(data = Pos_Sign,
               aes(x = pos, xend = pos,
                   y = -0.009, yend = 0.009,
                   color = "Diagnostic positions"),
               linewidth = 0.8) +
  
  # geom_text(aes(x = center, y = ifelse(type == "Exon",0.03,-0.03), label = paste(type,num)), size = 3,fontface = "bold",inherit.aes = FALSE)+

  scale_color_manual(values = c(
    "Exon" = "green3",
    "Intron" = "grey",
    "Diagnostic positions" = "red"
  )) +

  coord_cartesian(xlim = Region_Sign, ylim = c(-0.2, 0.2)) +
  scale_x_reverse(breaks = c(30558000, 30560000, 30562000)) +

  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=16),
        axis.title = element_text(size=14),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14)) +  

  labs(x = "Position",
       color = "Type")


ggsave(
       paste0(folder_out,"Gene_EDNRB_RegSign.png"),
       plot = p2,
       width = 10,
       height = 5,
       dpi = 150   
   )
