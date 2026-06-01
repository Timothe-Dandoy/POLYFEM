library(ggplot2)
library(dplyr)

file_fst    <- commandArgs(trailingOnly = TRUE)[1]
Diag_region <- commandArgs(trailingOnly = TRUE)[2]
species     <- commandArgs(trailingOnly = TRUE)[3]
prefix_plot <- commandArgs(trailingOnly = TRUE)[4]

L_col = list("P_guadarramae" = "darkblue","P_lusitanicus" = "gold3","P_bocagei"="green4","P_liolepis"="darkorange")
colour = L_col[[species]]

# Get data
df_FST <- read.table(
  file_fst,
  header = TRUE,
  sep = "",
  stringsAsFactors = FALSE
)

Diag_region = as.numeric(strsplit(Diag_region,"-")[[1]])

df_FST$fst    = as.numeric(df_FST$fst)
df_FST$midPos = as.numeric(df_FST$midPos)

#plot 

ymin = 0
ymax = 0.4
N_plots = 12

for (f in 1:N_plots){
  
  xmin = (f-1) * 5*10**6
  xmax =  f * 5*10**6
    
  p_fst_f = ggplot() +
    
    geom_rect(
      aes(xmin = Diag_region[1], xmax = Diag_region[2], ymin = -Inf, ymax = Inf, fill = "Diagnostic Region"),
      inherit.aes = FALSE,
      alpha = 0.5
    ) +
    
    geom_line(data = df_FST,aes(midPos, fst,color="Fst"), linewidth=0.5) +
    
    
    scale_color_manual(values = c("Fst" = colour)) +
    
    scale_fill_manual(values = c("Diagnostic Region" = "red")) +
    
    labs(title="Chromosome Z - FST",
         x="Physical Positions",
         y="Fst") +
    coord_cartesian(xlim = c(xmin, xmax),ylim = c(ymin,ymax))+
    
    theme_classic()
  
  
  ggsave(paste0(prefix_plot,f,'.png'), p_fst_f, width=10, height=6)
}


p_fst_all = ggplot() +
    
   geom_rect(
      aes(xmin = Diag_region[1], xmax = Diag_region[2], ymin = -Inf, ymax = Inf, fill = "Diagnostic Region"),
      inherit.aes = FALSE,
      alpha = 0.5
    ) +
    
    geom_line(data = df_FST,aes(midPos, fst,color="Fst"), linewidth=0.5) +
    
    
    scale_color_manual(values = c("Fst" = colour)) +
    
    scale_fill_manual(values = c("Diagnostic Region" = "red")) +
    
    labs(x="Positions",
         y="Fst") +
    
    coord_cartesian(ylim = c(ymin,ymax))+
    
    theme_classic()+
    theme(legend.position = "none")
  
  
ggsave(paste0(prefix_plot,'_all.png'), p_fst_all, width=70, height=5,limitsize = FALSE)


xmin = Diag_region[1]-100000
xmax =  Diag_region[2]+100000

p_fst_reg = ggplot() +
    
  geom_rect(
    aes(xmin = Diag_region[1], xmax = Diag_region[2], ymin = -Inf, ymax = Inf), 
    inherit.aes = FALSE,
    fill = NA,
    color = "red",
    alpha = 0.5
  ) +
    
  geom_step(data = df_FST,aes(midPos, fst), color = colour, linewidth=0.5) +
    
  labs(x="Positions",
       y="Fst") +
  
  coord_cartesian(xlim = c(xmin, xmax))+
    
  theme_classic()+
  theme(legend.position = "none")
  
  
ggsave(paste0(prefix_plot,'_region.png'), p_fst_reg, width=10, height=6)


# get the position with fst > 0.2

th = 0.2

df_Sign = df_FST[df_FST$fst > th,]
df_Sign = df_Sign[,c(3,5)]

write.table(
    df_Sign,
    file = paste0(prefix_plot,'_pos_sign.txt'),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)


