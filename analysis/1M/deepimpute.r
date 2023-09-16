# plot umap
library(umap)
library(rhdf5)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
set.seed(20)
dirs = "/home/suyanchi/project/dab/data/1M"
color = readRDS("/home/suyanchi/project/dab/color.rds")
color = color[c(-3,-4,-11)]

plot_umap <- function(df, methods_name, label){
  colourCount = length(unique(label))
  p = ggplot(df, aes(x, y, color = Cell)) +
    geom_point(size=1) +
    # ggtitle(methods_name) +
    # labs(x="UMAP_1", y="UMAP_2") +
    theme_bw() +
    scale_color_manual(values = colorRampPalette(color)(colourCount)) +
    theme(
      #axis.text.x=element_blank(),
      #axis.text.y=element_blank(),
      axis.ticks = element_blank(),
      # axis.title = element_blank(),
      axis.title.x = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(size = 0.5, colour = "black"),
      plot.title = element_text(hjust = 0.5))
  
  return(p)
}
load('/home/suyanchi/project/dab/results/1M/umap.RData')
# dca
h5file=H5Fopen('/home/suyanchi/project/dab/results/1M/deepimpute.h5')
dat = h5file$data
label = readRDS('/home/suyanchi/project/dab/data/1M/1M_cell.rds')
d.umap <- umap((dat))
df[[5]] <- data.frame(x = d.umap$layout[,1],
                           y = d.umap$layout[,2],
                           Cell = label)

p[[5]] = plot_umap(df, methods_name ='dropout', label)

save(df, p, file='/home/suyanchi/project/dab/results/1M/umap.RData')