library (ggplot2)

out_gwas   <- commandArgs(trailingOnly = TRUE)[1]
species    <- commandArgs(trailingOnly = TRUE)[2]
prefix_res <- commandArgs(trailingOnly = TRUE)[3]

Asso <- read.table(
  out_gwas,
  header = TRUE,
  sep = "",
  stringsAsFactors = FALSE
)

L_col = list("P_guadarramae" = "darkblue","P_lusitanicus" = "gold3","P_bocagei"="green4","P_liolepis"="darkorange")
colour = L_col[[species]]

Chromosome <- sub("NC_0413", "", Asso$Chromosome)
Chromosome <- sub("\\.1", "",Chromosome)
Asso$Num_chr <- as.numeric(Chromosome) -11

Asso <- Asso[!is.na(Asso$Num_chr),]
initial_Pos = Asso$Position
Seuil = 0.05/(1.5e9)
Seuil_log = -log10(Seuil)


for (i in 1:18){
  last_pos = max(Asso$Position[Asso$Num_chr==i],na.rm=TRUE)
  Asso$Position[Asso$Num_chr==i+1] = Asso$Position[Asso$Num_chr==i+1] +last_pos
}

Asso$Pval = pchisq(Asso$LRT, df=1, lower.tail=FALSE)
Asso$MinusLogPval = -log10(Asso$Pval)

p <- ggplot(
  Asso,
  aes(Position, MinusLogPval, color = factor(Num_chr %% 2))
) +
  geom_point(size = 0.3, alpha = 1) +
  geom_hline(yintercept = Seuil_log, color = "red") +
  scale_color_manual(values = c("grey40",colour)) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(size=16),
    axis.title = element_text(size=14),
    legend.title = element_text(size=16),
    legend.text = element_text(size=14)
  )+
  labs(y = expression(-log[10](pvalue)))+
  coord_cartesian(ylim = c(0, 50))+
  ggtitle(species)


 ggsave(
   paste0(prefix_res,"plot.png"),
   plot = p,
   width = 10,
   height = 5,
   dpi = 150   
 )

Asso_Sign = Asso[Asso$MinusLogPval > Seuil_log,]
Pos_sign = initial_Pos[Asso$MinusLogPval > Seuil_log]
Chr_sign = Asso$Num_chr[Asso$MinusLogPval > Seuil_log]
MLP_sign = Asso$MinusLogPval[Asso$MinusLogPval > Seuil_log]

sink(paste0(prefix_res,"diagnostic_positions.txt"))

for (i in 1:19) {
  
  pos <- Pos_sign[Chr_sign == i]
  mlp <- MLP_sign[Chr_sign == i]
  
  cat("\n============================\n")
  cat("Chromosome :", i, "\n")
  cat("Diagnostic Positions (Pval < 1e-9):\n")
  
  if (length(pos) == 0) {
    cat("  None\n")
  } else {
    cat(" ", paste(pos, collapse = ", "), "\n")
    cat("P values :\n")
    cat(" ", paste(mlp, collapse = ", "), "\n")
  }
}

sink()  


















