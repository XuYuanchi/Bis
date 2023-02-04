# seurat clustering for downsampled datasets

library(aricode)
library(mclustcomp)
library(R.matlab)
library(Seurat)
library(edgeR)
set.seed(0814)
dirs = "D:/project/Demo_Bis/clustering/"
color = readRDS("D:/project/Demo_Bulk/color.rds")

data_name = c("FetalBrain", "SmallIntestine", "Kidney", "Liver", 
              "Spleen", "Placenta", "FetalLiver", "Lung")
methods_name <- c("Dropout", "Bis", "ALRA", "DCA", "DeepImpute", "MAGIC", "SAVER", "scImpute", "SCRABBLE", "scScope", "scVI")
cal_metric <- function(count, label, gene){
  colnames(count) = label
  rownames(count) = c(1:nrow(count))
  x.seurat <- CreateSeuratObject(count)
  x.seurat <- NormalizeData(x.seurat)
  x.seurat <- FindVariableFeatures(x.seurat, selection.method = "vst", nfeatures = length(gene), verbose = FALSE)
  x.seurat <- ScaleData(x.seurat)
  
  x.seurat <- RunPCA(x.seurat, features = VariableFeatures(object = x.seurat))  
  x.seurat <- JackStraw(x.seurat, num.replicate = 100)
  x.seurat <- ScoreJackStraw(x.seurat, dims = 1:20)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
  x.seurat <- FindNeighbors(x.seurat, dims = 1:10)
  x.seurat <- FindClusters(x.seurat, resolution = 0.5)
  # 1:0.7
  # 2:0.3
  # 3:
  
  #ARI <- adjustedRandIndex(as.factor(label), x.seurat$seurat_clusters)
  ari <- ARI(as.factor(label), x.seurat$seurat_clusters)
  nmi <- NMI(as.factor(label), x.seurat$seurat_clusters)
  jaccard <- mclustcomp(as.integer(factor(label)), as.integer(x.seurat$seurat_clusters), types = "jaccard")$scores
  metric <- c(ari, nmi, jaccard)
  return (metric)
}
ARI = matrix(nrow = 11, ncol = 10)
NMI = matrix(nrow = 11, ncol = 10)
jaccard = matrix(nrow = 11, ncol = 10)
rownames(ARI) = methods_name
#colnames(ARI) = data_name
rownames(NMI) = methods_name
#colnames(NMI) = data_name
rownames(jaccard) = methods_name
#colnames(jaccard) = data_name

for (id in c(5)) {
  
  samp.file = paste0(dirs, 'data/', data_name[id], '_dropout.rds')
  Bis.file = list.files(paste0(dirs, 'Results/Bis/', data_name[id]), pattern = '0.9.mat', full.names = TRUE, recursive = TRUE)
  ALRA.file = list.files(paste0(dirs, 'Results/ALRA/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  DCA.file = list.files(paste0(dirs, 'Results/DCA/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  DeepImpute.file = list.files(paste0(dirs, 'Results/DeepImpute/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  MAGIC.file = list.files(paste0(dirs, 'Results/MAGIC/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  SAVER.file = list.files(paste0(dirs, 'Results/SAVER/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  scImpute.file = list.files(paste0(dirs, 'Results/scImpute/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  SCRABBLE.file = list.files(paste0(dirs, 'Results/SCRABBLE/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  scScope.file = list.files(paste0(dirs, 'Results/scScope/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  scVI.file = list.files(paste0(dirs, 'Results/scVI/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  
  
  
  dat <- vector("list", 11)
  
  dat[[1]] = as.matrix(readRDS(samp.file))
  cell_infor = readRDS(paste0(dirs, "data/", data_name[id], "_celltype.rds"))
  true_label <- as.vector(cell_infor$ClusterID)
  gene = rownames(dat[[1]])
  
  for (i in c(1:5)){
    dat[[2]] = readMat(Bis.file[i])[[1]]
    dat[[3]] = readRDS(ALRA.file[i])
    dat[[4]] = as.matrix(readRDS(DCA.file[i]))
    dat[[5]] = as.matrix(readRDS(DeepImpute.file[i]))
    dat[[6]] = as.matrix(readRDS(MAGIC.file[i]))
    dat[[7]] = t(as.matrix(readRDS(SAVER.file[i])))
    dat[[8]] = as.matrix(readRDS(scImpute.file[i]))
    # dat[[9]] = as.matrix(readMat(SCRABBLE.file[i])[[1]])
    dat[[9]] = as.matrix(readRDS(SCRABBLE.file))
    dat[[10]] = t(readRDS(scScope.file[i]))
    dat[[11]] = as.matrix(readRDS(scVI.file[i]))
    
    for (j in c(1:11)){
      
      metric = cal_metric(dat[[j]], true_label, gene)
      metric
      ARI[j, i] = metric[1]
      NMI[j, i] = metric[2]
      jaccard[j, i] = metric[3]
    }
    
  }
  
  save(file=paste0(dirs, 'Results/Bis/', data_name[id], '/metric_clustering_seurat.RData'), ARI=ARI, NMI=NMI, jaccard=jaccard)
  
}