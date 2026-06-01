library(ggplot2)
library("Biostrings")
library(patchwork)

Prefix_exon    <- commandArgs(trailingOnly = TRUE)[1]
Prefix_intron  <- commandArgs(trailingOnly = TRUE)[2]
Suffix         <- commandArgs(trailingOnly = TRUE)[3]
in_pheno       <- commandArgs(trailingOnly = TRUE)[4]
folder_out     <- commandArgs(trailingOnly = TRUE)[5]

# 1) Link Phenotype - Individuals
Phenotypes <- read.table(
  in_pheno,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)

colnames(Phenotypes) <- c("ID", "Phenotype")

IDs = Phenotypes$ID

Phenotypes$Phenotype[Phenotypes$ID == "T14419"] = NA           # Remove of the three individuals with missmatch phenotype in the haplotype network
Phenotypes$Phenotype[Phenotypes$ID == "T14457"] = NA
Phenotypes$Phenotype[Phenotypes$ID == "T14432"] = NA

Ligne = Phenotypes[Phenotypes$Phenotype == "S" & !is.na(Phenotypes$Phenotype),]
Male  = Phenotypes[Phenotypes$Phenotype == "RB" & !is.na(Phenotypes$Phenotype),]

# 2) Link Phenotype - Haplotype for all exons and introns 
List_out          = list()
List_out_diff     = list()
List_out_megadiff = list()


for (i in 1:15){
    # i even => intron  |  i odd => exon 
    Fasta = ifelse(i %% 2 == 1,paste0(Prefix_exon,i %/% 2 + 1,Suffix),paste0(Prefix_intron,i %/% 2,Suffix))

    seqs <- readDNAStringSet(Fasta)

    Sequences <- as.character(seqs)

    seq_list  <- strsplit(Sequences, split = "")
    names(seq_list) <- sub(":.*", "", names(seq_list))

    seq_lengths <- sapply(seq_list, length)
    max_length  <- max(seq_lengths)
    
    Maj_L <- c()
    Maj_M <- c()

    All_L <- c()
    All_M <- c()

    Nb_A_L   <- c()
    Nb_C_L   <- c()
    Nb_G_L   <- c()
    Nb_T_L   <- c()
    Nb_N_L   <- c()
    Nb_del_L <- c()

    Nb_A_M   <- c()
    Nb_C_M   <- c()
    Nb_G_M   <- c()
    Nb_T_M   <- c()
    Nb_N_M   <- c()
    Nb_del_M <- c()

    for (pos in 1:max_length){
        Nucs_L = c()
        for (el in Ligne$ID){
            seq = seq_list[[el]]
            if (length(seq) < pos  ){
                Nucs_L = c(Nucs_L,NA)
            } else{
                Nuc = seq[pos]
                Nucs_L = c(Nucs_L,Nuc)
            }
        }
        Nucs_M = c()
        for (el in Male$ID){
            seq = seq_list[[el]]
            if (length(seq) < pos  ){
                Nucs_M = c(Nucs_M,NA)
            } else{
                Nuc = seq[pos]
                Nucs_M = c(Nucs_M,Nuc)
            }
        }
        
        # Get the predominant nucleotide at the position, for each morph
        Maj_L[pos]    <- names(sort(table(Nucs_L), decreasing = TRUE))[1]
        Maj_M[pos]    <- names(sort(table(Nucs_M), decreasing = TRUE))[1]
    
        # 1 if all individuals of the same morphs have the same nucleotide at that position, 0 otherwise
        #All_L[pos] <- ifelse(length(table(Nucs_L[!(Nucs_L %in% c("N"))])) == 1,1,0)
        #All_M[pos] <- ifelse(length(table(Nucs_M[!(Nucs_M %in% c("N"))])) == 1,1,0)
        All_L[pos]    <- ifelse(length(table(Nucs_L[!is.na(Nucs_L)])) == 1,1,0)
        All_M[pos]    <- ifelse(length(table(Nucs_M[!is.na(Nucs_M)])) == 1,1,0)

        Nb_A_L[pos]   <- ifelse(is.na(as.numeric(table(Nucs_L)["A"])),0,as.numeric(table(Nucs_L)["A"]))
        Nb_C_L[pos]   <- ifelse(is.na(as.numeric(table(Nucs_L)["C"])),0,as.numeric(table(Nucs_L)["C"]))
        Nb_G_L[pos]   <- ifelse(is.na(as.numeric(table(Nucs_L)["G"])),0,as.numeric(table(Nucs_L)["G"]))
        Nb_T_L[pos]   <- ifelse(is.na(as.numeric(table(Nucs_L)["T"])),0,as.numeric(table(Nucs_L)["T"]))
        Nb_N_L[pos]   <- ifelse(is.na(as.numeric(table(Nucs_L)["N"])),0,as.numeric(table(Nucs_L)["N"]))
        Nb_del_L[pos] <- ifelse(is.na(as.numeric(table(Nucs_L)["-"])),0,as.numeric(table(Nucs_L)["-"]))

        Nb_A_M[pos]   <- ifelse(is.na(as.numeric(table(Nucs_M)["A"])),0,as.numeric(table(Nucs_M)["A"]))
        Nb_C_M[pos]   <- ifelse(is.na(as.numeric(table(Nucs_M)["C"])),0,as.numeric(table(Nucs_M)["C"]))
        Nb_G_M[pos]   <- ifelse(is.na(as.numeric(table(Nucs_M)["G"])),0,as.numeric(table(Nucs_M)["G"]))
        Nb_T_M[pos]   <- ifelse(is.na(as.numeric(table(Nucs_M)["T"])),0,as.numeric(table(Nucs_M)["T"]))
        Nb_N_M[pos]   <- ifelse(is.na(as.numeric(table(Nucs_M)["N"])),0,as.numeric(table(Nucs_M)["N"]))
        Nb_del_M[pos] <- ifelse(is.na(as.numeric(table(Nucs_M)["-"])),0,as.numeric(table(Nucs_M)["-"]))
    }
    
    out = ifelse(i %% 2 == 1,paste0("ex",i %/% 2 + 1),paste0("int",i %/% 2))

    OutPut <- data.frame(
    Position = 1:max_length,
    Ex_Int   = rep(out, max_length),
    Maj_L    = Maj_L,
    Maj_M    = Maj_M,
    All_L    = All_L,
    All_M    = All_M,
    Nb_A     = paste(Nb_A_L,"-",Nb_A_M),
    Nb_C     = paste(Nb_C_L,"-",Nb_C_M),
    Nb_G     = paste(Nb_G_L,"-",Nb_G_M),
    Nb_T     = paste(Nb_T_L,"-",Nb_T_M),
    Nb_N     = paste(Nb_N_L,"-",Nb_N_M),
    Nb_del   = paste(Nb_del_L,"-",Nb_del_M)
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
OutPut_diff     <- OutPut[OutPut$Maj_L != OutPut$Maj_M,]                               # When the majoritary alleles of ligne and male are different
OutPut_megadiff <- OutPut_diff[(OutPut_diff$All_L == 1 & OutPut_diff$All_M == 1),]     # When the majoritary alleles of ligne and male are different & all lignes have the same nucleotide ,and same for males


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

# 5) Genome Representation

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
       paste0(folder_out,"_Gene_EDNRB.png"),
       plot = p1 ,
       width = 10,
       height = 5,
       dpi = 150   
   )

ggsave(
       paste0(folder_out,"Gene_EDNRB_RegSign.png"),
       plot = p2,
       width = 10,
       height = 5,
       dpi = 150   
   )
