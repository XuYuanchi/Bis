# clustering 1M datasets
library(Seurat)
library(rhdf5)
library(aricode)
library(mclustcomp)
library(future)
options(future.globals.maxSize= 4000 * 1024^2)
plan("multiprocess", workers = 50)
dirs = '/home/suyanchi/project/dab/results/1M/'

cal_metric <- function(count, label){
  colnames(count) = c(1:ncol(count))
  rownames(count) = c(1:nrow(count))
  x.seurat <- CreateSeuratObject(count)
  x.seurat <- NormalizeData(x.seurat)
  x.seurat <- FindVariableFeatures(x.seurat, selection.method = "vst", nfeatures = nrow(count), verbose = FALSE)
  x.seurat <- ScaleData(x.seurat)
  
  x.seurat <- RunPCA(x.seurat, features = VariableFeatures(object=x.seurat))  
  # x.seurat <- JackStraw(x.seurat, num.replicate = 100)
  # x.seurat <- ScoreJackStraw(x.seurat, dims = 1:20)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
  x.seurat <- FindNeighbors(x.seurat, dims = 1:30)
  x.seurat <- FindClusters(x.seurat, resolution = 0.5)
  # 1:0.7
  # 2:0.3
  # 3:
  metric = vector('list', 4)
  x.seurat <- RunUMAP(x.seurat, dims = 1:30)
  metric[[4]] <- x.seurat@reductions$umap@cell.embeddings
  #ARI <- adjustedRandIndex(as.factor(label), x.seurat$seurat_clusters)
  metric[[1]] <- ARI(as.factor(label), x.seurat$seurat_clusters)
  metric[[2]] <- NMI(as.factor(label), x.seurat$seurat_clusters)
  metric[[3]] <- mclustcomp(as.integer(factor(label)), as.integer(x.seurat$seurat_clusters), types = "jaccard")$scores
  # metric <- c(ari, nmi, jaccard, x.embeding)
  return (metric)
}

# load label
label = readRDS('/home/suyanchi/project/dab/data/1M/1M_cell.rds')
# gene = readRDS('D:/project/dab/data/1M/1M_gene.rds')
ARI = matrix(nrow = 1, ncol = 6)
NMI = matrix(nrow = 1, ncol = 6)
jaccard = matrix(nrow = 1, ncol = 6)

h5file=H5Fopen('/home/suyanchi/project/dab/results/1M/alra.h5')
dat = h5file$data
metric = cal_metric(as.matrix(t(dat)), label)
saveRDS(file=paste0(dirs, 'metric_alra.rds'), metric)
h5closeAll()
print(metric[1:3])
print('alra done')