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
Diag_region   <- commandArgs(trailingOnly = TRUE)[3]
species       <- commandArgs(trailingOnly = TRUE)[4]
prefix_output <- commandArgs(trailingOnly = TRUE)[5]

L_col = list("P_guadarramae" = "darkblue","P_lusitanicus" = "gold3","P_bocagei"="green4","P_liolepis"="darkorange")
colour = L_col[[species]]

Phenotypes <- read.table(
  Pheno_ind,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)

colnames(Phenotypes) <- c("ID", "Phenotype")

#-----------------------------------------------------------------
# Scripts
#-----------------------------------------------------------------

striped = Phenotypes$ID[!is.na(Phenotypes$Phenotype) & Phenotypes$Phenotype == "S"] 
reticulated = Phenotypes$ID[!is.na(Phenotypes$Phenotype) & Phenotypes$Phenotype == "RB"]

# Get nucleotide sequences
ChrZ <- read.dna(input_fasta, format = "fasta", as.character=TRUE)
Diag_region = as.numeric(strsplit(Diag_region,"-")[[1]])

# 1) Compute genetic distance at each position between striped and reticulated blotched females
Nb_Pos <- ncol(ChrZ)
Dxy = numeric()

# We look at every positions
for (p in 1:Nb_Pos){
  p_striped = ChrZ[striped,p]
  p_reticulated = ChrZ[reticulated,p]

  p_striped = p_striped[p_striped != "n" & p_striped != "-"]
  p_reticulated = p_reticulated[p_reticulated != "n" & p_reticulated != "-"]

  #If one morph is not represente dat the position the distance is NA
  if (length(p_striped) == 0 || length(p_reticulated) == 0) {
    Dxy[p] <- NA
    next
  }
  
  # Compute the number of differences between positions of both morphs
  diff_matrix <- outer(p_striped, p_reticulated, FUN = "!=")

  N_diff = sum(diff_matrix)
  N_paires = length(diff_matrix)

  Dxy[p] = N_diff / N_paires
}

# 2) Average on 5000 bp windows

df_tmp = data.frame(pos = 1:Nb_Pos, Dxy = Dxy, window =0)

windows = (0:(Nb_Pos %/% 5000 + 1)) * 5000

for (i in 1:(length(windows)-1)){
  df_tmp$window[between(df_tmp$pos,windows[i],windows[i+1])] = i
}

Mean_windows = tapply(df_tmp$Dxy, df_tmp$window, mean, na.rm = TRUE)

idx_window = match(1:(length(windows)-1), unique(df_tmp$window))
Mean_all_windows = ifelse(is.na(idx_window), 0, Mean_windows[idx_window])

Dxy_window = data.frame(
  window = windows[2:length(windows)],
  Mean_Dxy = Mean_all_windows
)

# 3) plots 

# plot over 5000000 bp
N_plots = 12
for (f in 1:N_plots){
  
  xmin = (f-1) * 5*10**6
  xmax =  f * 5*10**6
  
  p_Dxy_f = ggplot() +
    
    geom_rect(
      aes(xmin = Diag_region[1], xmax = Diag_region[2], ymin = -Inf, ymax = Inf, fill = "Diagnostic Region"),
      inherit.aes = FALSE,
      alpha = 0.5
    ) +
    
    geom_step(data = Dxy_window,aes(window, Mean_Dxy,color="Dxy"), linewidth=0.5) +
    
    
    scale_color_manual(values = c("Dxy" = colour)) +
    
    scale_fill_manual(values = c("Diagnostic Region" = "red")) +
    
    labs(title="Chromosome Z - Mean Dxy over 5000 bp",
         x="Physical Positions",
         y="Dxy") +
    coord_cartesian(xlim = c(xmin, xmax))+
    
    theme_classic()
  
  ggsave(paste0(prefix_output,f,'.png'), p_Dxy_f, width=10, height=6)
}

# plot along the whole chromosome Z
p_Dxy_all = ggplot() +
  
  geom_rect(
    aes(xmin = Diag_region[1], xmax = Diag_region[2], ymin = -Inf, ymax = Inf, fill = "Diagnostic Region"),
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  
  geom_step(data = Dxy_window,aes(window, Mean_Dxy,color="Dxy"), linewidth=0.5) +
  
  
  scale_color_manual(values = c("Dxy" = colour)) +
  
  scale_fill_manual(values = c("Diagnostic Region" = "red")) +
  
  labs(title="Chromosome Z - Mean Dxy over 5000 bp",
       x="Physical Positions",
       y="Dxy") +
  
  theme_classic()

ggsave(paste0(prefix_output,'all.png'), p_Dxy_all, width=10, height=6)

# plot around the EDNRB region
xmin = Diag_region[1]-100000
xmax =  Diag_region[2]+100000

p_Dxy_reg = ggplot() +
  
  geom_rect(
    aes(xmin = Diag_region[1], xmax = Diag_region[2], ymin = -Inf, ymax = Inf, fill = "Diagnostic Region"),
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  
  geom_step(data = Dxy_window,aes(window, Mean_Dxy,color="Dxy"), linewidth=0.5) +
  
  
  scale_color_manual(values = c("Dxy" = colour)) +
  
  scale_fill_manual(values = c("Diagnostic Region" = "red")) +
  
  labs(title="Chromosome Z - Mean Dxy over 5000 bp",
       x="Physical Positions",
       y="Dxy") +
  
  coord_cartesian(xlim = c(xmin, xmax))+
  
  theme_classic()

ggsave(paste0(prefix_output,'region.png'), p_Dxy_reg, width=10, height=6)






