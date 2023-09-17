# plot umap with labels of different methods
set.seed(0814)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(aricode)
library(mclustcomp)
library(dplyr)
library(umap)

color = c(
  "#FED439E5", "#709AE1E5", "#8A9197E5", "#D2AF81E5", "#FD7446E5", "#D5E4A2E5", "#197EC0E5", "#F05C3BE5",
  "#46732EE5", "#71D0F5E5", "#370335E5"
)

cal_metric <- function(count, label, n_gene){
  colnames(count) = c(1:length(label))
  rownames(count) = c(1:nrow(count))
  x.seurat <- CreateSeuratObject(count)
  x.seurat <- NormalizeData(x.seurat)
  x.seurat <- FindVariableFeatures(x.seurat, selection.method = "vst", nfeatures = n_gene, verbose = FALSE)
  x.seurat <- ScaleData(x.seurat)
  
  x.seurat <- RunPCA(x.seurat, features = VariableFeatures(object = x.seurat))  
  x.seurat <- JackStraw(x.seurat, num.replicate = 100)
  x.seurat <- ScoreJackStraw(x.seurat, dims = 1:20)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
  x.seurat <- FindNeighbors(x.seurat, dims = 1:10)
  x.seurat <- FindClusters(x.seurat, resolution = 0.5)
  # x.seurat <- RunUMAP(x.seurat)
  # 1:0.7
  # 2:0.3
  # 3:
  
  #ARI <- adjustedRandIndex(as.factor(label), x.seurat$seurat_clusters)
  ari <- ARI(as.factor(label), x.seurat$seurat_clusters)
  nmi <- NMI(as.factor(label), x.seurat$seurat_clusters)
  jaccard <- mclustcomp(as.integer(factor(label)), as.integer(x.seurat$seurat_clusters), types = "jaccard")$scores
  metric <- c(ari, nmi, jaccard)
  # umap = x.seurat@reductions$umap@cell.embeddings %>%
    #as.data.frame() %>% 
    #cbind(tx = x.seurat$seurat_clusters) # tx可替换为你所希望的列名
  results = list(metric, x.seurat$seurat_clusters)
  
  return (results)
}

plot_umap <- function(df, methods, color){
  
  p=ggplot(df, aes_string("UMAP_1", "UMAP_2", color=methods))+
    geom_point(size = 0.05) +
    ggtitle(methods) +
    # labs(x="UMAP_1", y="UMAP_2") +
    theme_bw() +
    scale_color_manual(values = color) +
    # scale_color_manual(values = color1) +
    theme(
      # axis.text.x=element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      # axis.title.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      # panel.border = element_blank(),
      axis.line = element_line(size = 0.5, colour = "black"),
      legend.position = "bottom",
      # legend.key.size=unit(25,'pt')
      # plot.title = element_text(hjust = 0.5)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3) ) ) +
    geom_segment(aes(x = min(df$UMAP_1) , y = min(df$UMAP_2) ,
                     xend = min(df$UMAP_1)+3, yend = min(df$UMAP_2) ),
                 colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm")))+ 
    geom_segment(aes(x = min(df$UMAP_1)  , y = min(df$UMAP_2)  ,
                     xend = min(df$UMAP_1) , yend = min(df$UMAP_2) + 3),
                 colour = "black", size=0.5,arrow = arrow(length = unit(0.3,"cm"))) +
    annotate("text", x = min(df$UMAP_1) +1, y = min(df$UMAP_2) -1, label = "UMAP_1",
             color="black",size = 2 ) + 
    annotate("text", x = min(df$UMAP_1) -1, y = min(df$UMAP_2) + 1, label = "UMAP_2",
             color="black",size = 2, angle=90) 
  return(p)
}


# load data
dirs = "D:/project/Demo_Bis/clustering/"
color = readRDS("D:/project/Demo_Bulk/color.rds")

data_name = c("FetalBrain", "SmallIntestine", "Kidney", "Liver", 
              "Spleen", "Placenta", "FetalLiver", "Lung")
methods_name <-  c("True", "Dropout", "Bis", "ALRA", "DCA", "DeepImpute", "MAGIC", "SAVER", "scImpute", "SCRABBLE", "scScope", "scVI")

for(i in c(8)){
  
  cell_infor = readRDS(paste0(dirs, "data/", data_name[i], "_celltype.rds"))
  
  true_label <- as.vector(cell_infor$Annotation)
  true_label = gsub("\\(.*?\\)", "", true_label)

  samp.file = paste0(dirs, 'data/', data_name[i], '_dropout.rds')

  dat = as.matrix(readRDS(samp.file))
  
  n_gene = dim(dat)[1]
  
  # re = cal_metric(dat, true_label, n_gene)
  
  # drop_label = re[[2]]
  
  d.umap <- umap(t(dat))
  dat_umap <- data.frame(UMAP_1 = d.umap$layout[,1],
                             UMAP_2 = d.umap$layout[,2],
                         # Dropout = drop_label,
                         True_Label=true_label)

  
  # load labels of different methods
  load(paste0(dirs, 'Results/Bis/', data_name[i], '/labels_clustering_seurat_later.RData'))
  
  dat_umap = cbind(dat_umap, pre_label)
  
  colnames(dat_umap) = c("UMAP_1", "UMAP_2","True", "Dropout", "Bis", "ALRA", "DCA", "DeepImpute", "MAGIC", "SAVER", "scImpute", "SCRABBLE", "scScope", "scVI")
  p=list()
  colourCount = length(unique(dat_umap$True))
  color1 = colorRampPalette(color)(colourCount)
  for(j in c(1:12)){
    p[[j]]=plot_umap(dat_umap, methods_name[j],color1)
    
  }
  saveRDS(p, paste0(dirs,  'plots/', data_name[i], '_anno.rds'))
  
}

for(i in c(8)){
  p=readRDS(paste0(dirs,  'plots/', data_name[i], '_anno.rds'))
  for(j in c(1:12)){
    p[[j]] = p[[j]] + theme(legend.text = element_text(size = 3),
                            legend.title = element_blank(),
                            legend.position = "none"
                            )
      # guides(color=guide_legend(override.aes = list(size=1)))
    ggsave(p[[j]], filename=paste0(dirs,  'plots/anno_umap/', data_name[i],'_',methods_name[j], '_anno.pdf'),
           width = 5, height = 5)
  }
}
# 3, 4, 5
# 4,5,7
# 3,4,5
# 3,5,7
