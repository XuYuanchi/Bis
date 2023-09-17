# plot barplot
dirs = "D:/project/Demo_Bis/clustering/Results/Bis/"
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
data_name = c("Spleen", "Kidney", "FetalLiver", "Lung")
color = readRDS("D:/project/Demo_Bulk/color.rds")
df = data.frame()
for (i in c(1:4)) {
  load(paste0(dirs, data_name[i], '/', 'metric_clustering_seurat.RData'))
  # delete NA
  ARI = data.frame(t(ARI[,1:5]))
  NMI = data.frame(t(NMI[,1:5]))
  jaccard = data.frame(t(jaccard[,1:5]))
  #if(i==4){
  #  ARI = ARI[,-2]
  #  NMI = NMI[,-2]
  #  jaccard = jaccard[,-2]}
  
  ARI$metric = rep('ARI', times=5)
  NMI$metric = rep('NMI', times=5)
  jaccard$metric = rep('jaccard', times=5)
  
  ARI$dataname = rep(data_name[i], times=5)
  NMI$dataname = rep(data_name[i], times=5)
  jaccard$dataname = rep(data_name[i], times=5)
  
  
  ARI = gather(ARI, Methods, value, -metric, -dataname)
  NMI = gather(NMI, Methods, value, -metric, -dataname)
  jaccard = gather(jaccard, Methods, value, -metric, -dataname)
  
  df = rbind(df, ARI)
  df = rbind(df, NMI)
  df = rbind(df, jaccard)
}
df$Methods = factor(df$Methods, levels = c("Dropout", "Bis", "ALRA", "DCA", "DeepImpute", "MAGIC", "SAVER", 
                                           "scImpute","SCRABBLE", "scScope", "scVI"))
df$metric = factor(df$metric, levels = c("ARI", "NMI","jaccard"))

df$dataname = factor(df$dataname, levels = c("Spleen", "Kidney", "FetalLiver", "Lung"))

p = ggbarplot(df, x = "Methods", y = "value", fill = "Methods", palette = color, 
              width = 0.8, ylab = "value", position = position_dodge(0.8),
              add = c('mean_se'),facet.by = c("metric","dataname"), 
              panel.labs = list(dataname=c("Spleen", "Placenta", "FetalLiver", "Lung")), ggtheme = theme_bw()) + rotate_x_text(45)

ggpar(p, ylim = c(0.2,0.85))
saveRDS(p, file = 'D:/project/Demo_Bis/clustering/plots/metric.rds')
ggsave(plot = p, filename = 'D:/project/Demo_Bis/clustering/plots/metric.pdf', height = 8, width = 14)

  
p = facet(p, facet.by = c("metric","dataname"),
          scales = "free_y",
          panel.labs = list(
            metric=c("ARI", "NMI","jaccard"),
            dataname=c('Spleen', 'Kidney', 'FetalLiver', 'Lung')),
          panel.labs.font = list(
            size = 14,
            face = "bold"))
