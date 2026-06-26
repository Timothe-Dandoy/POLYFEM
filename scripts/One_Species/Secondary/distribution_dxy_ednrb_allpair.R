#-----------------------------------------------------------------
# Packages
#-----------------------------------------------------------------

library(ape)
library(ggplot2)
library(dplyr)

#-----------------------------------------------------------------
# Data
#-----------------------------------------------------------------

input_fasta   <- commandArgs(trailingOnly = TRUE)[1]
Pheno_ind     <- commandArgs(trailingOnly = TRUE)[2]
species       <- commandArgs(trailingOnly = TRUE)[3]
prefix_output <- commandArgs(trailingOnly = TRUE)[4]

L_col = list("P_guadarramae" = "darkblue","P_lusitanicus" = "gold3","P_bocagei"="green4","P_liolepis"="darkorange")
colour = L_col[[species]]

#-----------------------------------------------------------------
# Scripts
#-----------------------------------------------------------------

# Get nucleotide sequence of the EDNRB gene
ChrZ <- read.dna(input_fasta, format = "fasta", as.character=TRUE)

# 1) Compute genetic distance between each pair of individuals
Nb_Pos <- ncol(ChrZ)
Nb_ind <- nrow(ChrZ)
Dxy = matrix(data = NA, nrow = Nb_ind, ncol = Nb_ind)

#look at all combinaison
for (i in 1:(Nb_ind-1)){
    for (j in (i+1):Nb_ind){
        diff = 0
        bp   = 0
        #compared the two individuals at each position
        for (p in 1:Nb_Pos){
            nuc_1 <- ChrZ[i,p]
            nuc_2 <- ChrZ[j,p]
            if ((nuc_1 != "n") & (nuc_1 != "-") & (nuc_2 != "n") & (nuc_2 != "-")){
                bp <- bp + 1
                if (nuc_1 != nuc_2){
                    diff <- diff + 1
                }
            }
        }
        Dxy[i,j] = diff / bp
        Dxy[j,i] = diff / bp
    }
}

diag(Dxy) = 0
uptri_Dxy <- Dxy[upper.tri(Dxy)]

# 2) Distribution histogramme of Dxy
png(filename = paste0(prefix_output,'_hist_ind.png'), width=1000, height=600)

hist(uptri_Dxy, col=colour, main="Dxy histogram - EDNRB", xlab="dXY", ylab="Count")

dev.off()

# 3) general table => Dxy value of each pair
Dxy_df = data.frame(Dxy)
rownames(Dxy_df) <- sub(":.*", "", rownames(ChrZ))
colnames(Dxy_df) <- sub(":.*", "", rownames(ChrZ))

write.table(
  Dxy_df,
  file = paste0(prefix_output,"_table.txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

# List of phenotypes
Phenotypes <- read.table(
  Pheno_ind,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)

colnames(Phenotypes) <- c("ID", "Phenotype")

Striped = Phenotypes$ID[!is.na(Phenotypes$Phenotype) & Phenotypes$Phenotype == "S"] 
Reticulated = Phenotypes$ID[!is.na(Phenotypes$Phenotype) & Phenotypes$Phenotype == "RB"]

# 4) Sub tables => Dxy value of each pair based on phenotype (S vs RB ; S vs S and RB vs RB)
df_Dxy_S_RB  = Dxy_df[Striped,Reticulated]
df_Dxy_S_S   = Dxy_df[Striped,Striped]
df_Dxy_RB_RB = Dxy_df[Reticulated,Reticulated]

write.table(
  df_Dxy_S_RB,
  file = paste0(prefix_output,"_table_S_RB.txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

write.table(
  df_Dxy_S_S,
  file = paste0(prefix_output,"_table_S_S.txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

write.table(
  df_Dxy_RB_RB,
  file = paste0(prefix_output,"table_RB_RB.txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

# 5) Compute the mean Dxy of couples between S and S, S and RB or RB and RB
Dxy_S_RB = as.matrix(df_Dxy_S_RB)
Dxy_S_S = as.matrix(df_Dxy_S_S)
Dxy_RB_RB = as.matrix(df_Dxy_RB_RB)

uptri_Dxy_S_S <- Dxy_l_l[upper.tri(Dxy_S_S)]
uptri_Dxy_RB_RB <- Dxy_t_t[upper.tri(Dxy_RB_RB)]

Mean_S_RB = mean(Dxy_S_RB, na.rm = TRUE)
Mean_S_S = mean(uptri_Dxy_S_S, na.rm = TRUE)
Mean_RB_RB = mean(uptri_Dxy_RB_RB, na.rm = TRUE)

sd_S_RB = sd(Dxy_S_RB, na.rm = TRUE)
sd_S_S = sd(uptri_Dxy_S_S, na.rm = TRUE)
sd_RB_RB = sd(uptri_Dxy_RB_RB, na.rm = TRUE)

stat_dxy = data.frame(
    Comparison = c("S x RB","S x S","RB x RB" ),
    Mean       = c(Mean_S_RB,Mean_S_S,Mean_RB_RB),
    SD         = c(sd_S_RB,sd_S_S,sd_RB_RB)
)

write.table(
  stat_dxy,
  file = paste0(prefix_output,"_table_mean.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# 6) Distribution histogramme of Dxy coloured by couple type

vec_Dxy_S_RB = as.vector(Dxy_S_RB)

len_S_RB = length(vec_Dxy_S_RB)
len_S_S = length(uptri_Dxy_S_S)
len_RB_RB = length(uptri_Dxy_RB_RB)

all_Dxy_sorted = c(vec_Dxy_S_RB,uptri_Dxy_S_S,uptri_Dxy_RB_RB)

couple_type = c(rep("S - RB",len_S_RB),rep("S - S",len_S_S),rep("RB - RB",len_RB_RB))

df_col_couple = data.frame(
  Dxy  = all_Dxy_sorted,
  couple_type = couple_type
)

g_dxy_ednrb_col <- ggplot(data = df_col_couple, aes(x = Dxy, fill = couple_type)) +
    
    geom_histogram(bins = 30, alpha = 1, position = "stack", color = "black")+

    scale_fill_manual(values = c("S - RB" = "#1f77b4",
                                 "S - S" = "#d62728",
                                 "RB - RB" = "#2ca02c")) +
    
    labs(title="Dxy distribution - EDNRB",
         x="Dxy",
         y="Count",
         fill = "Couple Type") +
    
    theme_classic()

ggsave(paste0(prefix_output,'_hist_ind_col_couple.png'), g_dxy_ednrb_col, width=10, height=6)


