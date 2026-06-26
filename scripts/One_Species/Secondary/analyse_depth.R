#-----------------------------------------------------------------
# packages
#-----------------------------------------------------------------

library(readxl)
library(tidyr)
library(ggplot2)

#-----------------------------------------------------------------
# Data
#-----------------------------------------------------------------

file_depth  <- commandArgs(trailingOnly = TRUE)[1]
bam_list    <- commandArgs(trailingOnly = TRUE)[2]
table_pheno <- commandArgs(trailingOnly = TRUE)[3]
prefix      <- commandArgs(trailingOnly = TRUE)[4]
diagnostic  <- commandArgs(trailingOnly = TRUE)[5] 
species    <- commandArgs(trailingOnly = TRUE)[6] 

# Read all table
bams <- readLines(bam_list)
IDs <- sub("\\.bam$", "", basename(bams))

Depth <- read.table(
  file_depth,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)

colnames(Depth) <- c("CHR", "Position",IDs)

BD = read_xlsx(path =table_pheno)
BD <- BD[BD$Species == species,]
BD <- BD[!is.na(BD$'genetic data') & BD$'genetic data' == "Y",]

BD$Phenotype[BD$Sample_ID == "T14419"] = NA          
BD$Phenotype[BD$Sample_ID == "T14457"] = NA
BD$Phenotype[BD$Sample_ID == "T14432"] = NA

IDs_S  = BD$Sample_ID[BD$Phenotype == "S"  & !is.na(BD$Phenotype)]
IDs_RB = BD$Sample_ID[BD$Phenotype == "RB" & !is.na(BD$Phenotype)]

lim_diagnostic = as.numeric(strsplit(diagnostic, "-")[[1]])

#-----------------------------------------------------------------
# Scripts
#-----------------------------------------------------------------

# Depth per phenotype

Depth$D_mean_S <- rowMeans(Depth[, IDs_S], na.rm = TRUE)
Depth$D_sd_S   <- apply(Depth[, IDs_S], 1, sd, na.rm = TRUE)

Depth$D_mean_RB <- rowMeans(Depth[, IDs_RB], na.rm = TRUE)
Depth$D_sd_RB   <- apply(Depth[, IDs_RB], 1, sd, na.rm = TRUE)


# Plots

p_S_RB <- ggplot(Depth, aes(x = Position)) +
  
  geom_line(aes(y = D_mean_S, color = "Depth_S"), size = 1) + 
  
  geom_line(aes(y = D_mean_RB, color = "Depth_RB"), size = 1) + 
  
  geom_vline(aes(xintercept = lim_diagnostic[1], color = "Hypothetic Invertion Limits"),linetype = "dashed",size = 0.5 )+
  
  geom_vline(aes(xintercept = lim_diagnostic[2], color = "Hypothetic Invertion Limits"),linetype = "dashed",size = 0.5 )+
  
  scale_color_manual(
    name = "Species",
    values = c("Depth_S" = "red3",
               "Depth_RB" = "green3",
               "Hypothetic Invertion Limits" = "blue")
  ) +
  
  labs(title = paste(species, "- S and RB"),
       x = "Genomic Position",
       y = "Depth") +

  scale_x_continuous(n.breaks = 10) +
  coord_cartesian(ylim = c(0, 50))+

  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

ggsave(file = paste0(prefix,"S_RB.png"), plot = p_S_RB, dpi = 150, width = 14, height = 6)


p_S <- ggplot(Depth, aes(x = Position,y = D_mean_S)) +
  geom_line(color = "red3", size = 1) + 
  geom_ribbon(aes(ymin = D_mean_S - D_sd_S, 
                  ymax = D_mean_S + D_sd_S),
              fill = "pink", alpha = 0.3) + 

  geom_vline(xintercept = lim_diagnostic[1],color = "blue",linetype = "dashed",size = 0.5 )+

  geom_vline(xintercept = lim_diagnostic[2],color = "blue",linetype = "dashed",size = 0.5 )+
  
  labs(title = paste(species, "- S"),
       x = "Genomic Position",
       y = "Depth") +

  scale_x_continuous(n.breaks = 10) +
  coord_cartesian(ylim = c(0, 50))+

  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

ggsave(file = paste0(prefix,"S.png"), plot = p_S, dpi = 150, width = 14, height = 6)


p_RB <- ggplot(Depth, aes(x = Position,y = D_mean_RB)) +
  geom_line(color = "green3", size = 1) + 
  geom_ribbon(aes(ymin = D_mean_RB - D_sd_RB, 
                  ymax = D_mean_RB + D_sd_RB),
              fill = "lightgreen", alpha = 0.3) + 

  geom_vline(xintercept = lim_diagnostic[1],color = "blue",linetype = "dashed",size = 0.5 )+

  geom_vline(xintercept = lim_diagnostic[2],color = "blue",linetype = "dashed",size = 0.5 )+
  
  labs(title = paste(species, "- RB"),
       x = "Genomic Position",
       y = "Depth") +

  scale_x_continuous(n.breaks = 10) +
  coord_cartesian(ylim = c(0, 50))+

  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

ggsave(file = paste0(prefix,"RB.png"), plot = p_RB, dpi = 150, width = 14, height = 6)