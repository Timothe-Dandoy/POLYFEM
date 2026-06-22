#-----------------------------------------------------------------
# Packages
#-----------------------------------------------------------------

library(ggplot2)
library(dplyr)

#-----------------------------------------------------------------
# Data
#-----------------------------------------------------------------

file_ld     <- commandArgs(trailingOnly = TRUE)[1]
file_bim    <- commandArgs(trailingOnly = TRUE)[2]
Diag_region <- commandArgs(trailingOnly = TRUE)[3]
window_mean <- commandArgs(trailingOnly = TRUE)[4]
species     <- commandArgs(trailingOnly = TRUE)[5]
prefix_plot <- commandArgs(trailingOnly = TRUE)[6]

L_col = list("P_guadarramae" = "darkblue","P_lusitanicus" = "gold3","P_bocagei"="green4","P_liolepis"="darkorange")
colour = L_col[[species]]

df_LD <- read.table(
  file_ld,
  header = TRUE,
  sep = "",
  stringsAsFactors = FALSE
)

BIM <- read.table(
  file_bim,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)

df_LD$BP_B = as.numeric(df_LD$BP_B)

All_Pos = as.numeric(BIM[,4])
Diag_region = as.numeric(strsplit(Diag_region,"-")[[1]])

#-----------------------------------------------------------------
# Scripts
#-----------------------------------------------------------------

# 1) get LD for all adjacent SNP

idx = match(All_Pos, df_LD$BP_B)
LD = ifelse(is.na(idx), 0, df_LD$R2[idx])

# 2) Mean LD over window_mean base pair

window_mean = as.numeric(window_mean)
df_tmp = data.frame(pos = All_Pos, LD = LD, window =0)

M_Pos = max(All_Pos)
windows = (0:(M_Pos %/% window_mean + 1)) * window_mean

for (i in 1:(length(windows)-1)){
  df_tmp$window[between(df_tmp$pos,windows[i],windows[i+1])] = i
}

Mean_windows = tapply(df_tmp$LD, df_tmp$window, mean)

idx_window = match(1:(length(windows)-1), unique(df_tmp$window))
Mean_all_windows = ifelse(is.na(idx_window), 0, Mean_windows[idx_window])

df_final = data.frame(
  window = windows[2:length(windows)],
  Mean_LD = Mean_all_windows
  )


# 3) plots 

N_plots = 12

# plot over 5000000 bp 
for (f in 1:N_plots){
  
  xmin = (f-1) * 5*10**6
  xmax =  f * 5*10**6
    
  p_LD_f = ggplot() +
    
    geom_rect(
      aes(xmin = Diag_region[1], xmax = Diag_region[2], ymin = -Inf, ymax = Inf, fill = "Diagnostic Region"),
      inherit.aes = FALSE,
      alpha = 1
    ) +
    
    geom_step(data = df_final,aes(window, Mean_LD,color="mean_R2"), linewidth=0.5) +
    
    
    scale_color_manual(values = c("mean_R2" = colour)) +
    
    scale_fill_manual(values = c("Diagnostic Region" = "red")) +
    
    labs(title="Chromosome Z - Mean LD over 5000 bp",
         x="Physical Positions",
         y="R2") +
    coord_cartesian(xlim = c(xmin, xmax))+
    
    theme_classic()+
    theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(size=16),
    axis.title = element_text(size=14),
    legend.title = element_text(size=16),
    legend.text = element_text(size=14)
  )
  
  
  ggsave(paste0(prefix_plot,f,'.png'), p_LD_f, width=10, height=6)
}

# plot along the whole chromosome Z   
p_LD_all = ggplot() +
    
   geom_rect(
      aes(xmin = Diag_region[1], xmax = Diag_region[2], ymin = -Inf, ymax = Inf, fill = "Diagnostic Region"),
      inherit.aes = FALSE,
      alpha = 1
    ) +
    
    geom_step(data = df_final,aes(window, Mean_LD,color="mean_R2"), linewidth=0.5) +
    
    
    scale_color_manual(values = c("mean_R2" = colour)) +
    
    scale_fill_manual(values = c("Diagnostic Region" = "red")) +
    
    labs(x="Physical Positions", y="R2") +
    
    theme_classic()+
    theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(size=16),
    axis.title = element_text(size=14),
    legend.title = element_text(size=16),
    legend.text = element_text(size=14)
  )
  
  
ggsave(paste0(prefix_plot,'_all.png'), p_LD_all, width=70, height=5,limitsize = FALSE)

# plot around the EDNRB region
xmin = Diag_region[1]-100000
xmax =  Diag_region[2]+100000

p_LD_reg = ggplot() +
    
  geom_rect(
    aes(xmin = Diag_region[1], xmax = Diag_region[2], ymin = -Inf, ymax = Inf), 
    inherit.aes = FALSE,
    fill = NA,
    color = "red",
    alpha = 0.5
  ) +
    
  geom_step(data = df_final,aes(window, Mean_LD), color = colour, linewidth=0.5) +
    
  labs(x="Physical Positions",
       y="R2") +
  
  coord_cartesian(xlim = c(xmin, xmax))+
     
  theme_classic()+
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(size=16),
    axis.title = element_text(size=14),
    legend.title = element_text(size=16),
    legend.text = element_text(size=14)
  )
  
ggsave(paste0(prefix_plot,'_region.png'), p_LD_reg, width=10, height=6)

