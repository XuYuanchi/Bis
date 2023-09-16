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
p = vector("list", 7)
df = vector('list', 7)
# dropout 
h5file=H5Fopen('/home/suyanchi/project/dab/data/1M/1M.h5')
dat = h5file$data
label = colnames(dat)
d.umap <- umap(t(dat))
df[[1]] <- data.frame(x = d.umap$layout[,1],
                           y = d.umap$layout[,2],
                           Cell = label)

p[[1]] = plot_umap(df, methods_name ='dropout', label)

# bis 
h5file=H5Fopen('/home/suyanchi/project/dab/results/1M/bis.h5')
dat = h5file$data

d.umap <- umap(t(dat))
df[[2]] <- data.frame(x = d.umap$layout[,1],
                 y = d.umap$layout[,2],
                 Cell = label)

p[[2]]=plot_umap(df, methods_name ='Bis', p = label)

# alra 
dat = readRDS('/home/suyanchi/project/dab/results/1M/alra.rds')

d.umap <- umap(t(dat))
df[[3]] <- data.frame(x = d.umap$layout[,1],
                 y = d.umap$layout[,2],
                 Cell = label)

p[[3]] = plot_umap(df, methods_name ='ALRA', label)

# dca 
h5file=H5Fopen('/home/suyanchi/project/dab/results/1M/dca.h5')
dat = h5file$data

d.umap <- umap((dat))
df[[4]] <- data.frame(x = d.umap$layout[,1],
                 y = d.umap$layout[,2],
                 Cell = label)

p[[4]] = plot_umap(df, methods_name ='DCA', label)

# deepimpute 
h5file=H5Fopen('/home/suyanchi/project/dab/results/1M/deepimpute.h5')
dat = h5file$data

d.umap <- umap((dat))
df[[5]] <- data.frame(x = d.umap$layout[,1],
                 y = d.umap$layout[,2],
                 Cell = label)

p[[5]] = plot_umap(df, methods_name ='DeepImpute', label)

# scscope 
h5file=H5Fopen('/home/suyanchi/project/dab/results/1M/scscope.h5')
dat = h5file$data

d.umap <- umap((dat))
df[[6]] <- data.frame(x = d.umap$layout[,1],
                 y = d.umap$layout[,2],
                 Cell = label)

p[[6]] = plot_umap(df, methods_name ='scScope', label)

# scvi 
h5file=H5Fopen('/home/suyanchi/project/dab/results/1M/scvi.h5')
dat = h5file$data

d.umap <- umap((dat))
df[[7]] <- data.frame(x = d.umap$layout[,1],
                 y = d.umap$layout[,2],
                 Cell = label)

p[[7]] = plot_umap(df, methods_name ='scVI', label)

save(p, df, file = '/home/suyanchi/project/dab/results/1M/umap.RData')
