library(ape)
library(ggplot2)
library(dplyr)

input_fasta   <- commandArgs(trailingOnly = TRUE)[1]
Pheno_ind     <- commandArgs(trailingOnly = TRUE)[2]
species       <- commandArgs(trailingOnly = TRUE)[3]
prefix_output <- commandArgs(trailingOnly = TRUE)[4]

L_col = list("P_guadarramae" = "darkblue","P_lusitanicus" = "gold3","P_bocagei"="green4","P_liolepis"="darkorange")
colour = L_col[[species]]

# Get nucleotide sequences
ChrZ <- read.dna(input_fasta, format = "fasta", as.character=TRUE)

# Compute genetic distance between each pair of individuals
Nb_Pos <- ncol(ChrZ)
Nb_ind <- nrow(ChrZ)
Dxy = matrix(data = NA, nrow = Nb_ind, ncol = Nb_ind)

for (i in 1:(Nb_ind-1)){
    for (j in (i+1):Nb_ind){
        diff = 0
        bp   = 0
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
# plot

uptri_Dxy <- Dxy[upper.tri(Dxy)]

png(filename = paste0(prefix_output,'_hist_ind.png'), width=1000, height=600)

hist(uptri_Dxy, col=colour, main="Dxy histogram - EDNRB", xlab="dXY", ylab="Count")

dev.off()

# general table

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

Lignees = Phenotypes$ID[!is.na(Phenotypes$Phenotype) & Phenotypes$Phenotype == "S"] 
Tachees = Phenotypes$ID[!is.na(Phenotypes$Phenotype) & Phenotypes$Phenotype == "RB"]

# Sub tables

df_Dxy_l_t = Dxy_df[Lignees,Tachees]
df_Dxy_l_l = Dxy_df[Lignees,Lignees]
df_Dxy_t_t = Dxy_df[Tachees,Tachees]

write.table(
  df_Dxy_l_t,
  file = paste0(prefix_output,"_table_S_RB.txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

write.table(
  df_Dxy_l_l,
  file = paste0(prefix_output,"_table_S_S.txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

write.table(
  df_Dxy_t_t,
  file = paste0(prefix_output,"table_RB_RB.txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

# Mean 
Dxy_l_t = as.matrix(df_Dxy_l_t)
Dxy_l_l = as.matrix(df_Dxy_l_l)
Dxy_t_t = as.matrix(df_Dxy_t_t)

uptri_Dxy_l_l <- Dxy_l_l[upper.tri(Dxy_l_l)]
uptri_Dxy_t_t <- Dxy_t_t[upper.tri(Dxy_t_t)]

Mean_l_t = mean(Dxy_l_t, na.rm = TRUE)
Mean_l_l = mean(uptri_Dxy_l_l, na.rm = TRUE)
Mean_t_t = mean(uptri_Dxy_t_t, na.rm = TRUE)

sd_l_t = sd(Dxy_l_t, na.rm = TRUE)
sd_l_l = sd(uptri_Dxy_l_l, na.rm = TRUE)
sd_t_t = sd(uptri_Dxy_t_t, na.rm = TRUE)

stat_dxy = data.frame(
    Comparison = c("S x RB","S x S","RB x RB" ),
    Mean       = c(Mean_l_t,Mean_l_l,Mean_t_t),
    SD         = c(sd_l_t,sd_l_l,sd_t_t)
)

write.table(
  stat_dxy,
  file = paste0(prefix_output,"_table_mean.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# Plot general 1 color by couple type

vec_Dxy_l_t = as.vector(Dxy_l_t)

len_l_t = length(vec_Dxy_l_t)
len_l_l = length(uptri_Dxy_l_l)
len_t_t = length(uptri_Dxy_t_t)

all_Dxy_sorted = c(vec_Dxy_l_t,uptri_Dxy_l_l,uptri_Dxy_t_t)

couple_type = c(rep("S - RB",len_l_t),rep("S - S",len_l_l),rep("RB - RB",len_t_t))

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


