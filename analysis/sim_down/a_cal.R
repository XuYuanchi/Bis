dirs = "D:/project/Demo_Bis/downsample/"
color = readRDS("D:/project/Demo_Bulk/color.rds")
# load library
library(R.matlab)
library(ggplot2)
library(edgeR)
source(paste0(dirs, 'scripts/a_cal_lib.R'))

data_name = c('Liver', 'Lung', 'Heart', 'Marrow')
# read filename
for(id in c(1)){
  
  samp.file = list.files(paste0(dirs, 'datasets/', data_name[id]), pattern = '._samp.rds', full.names = TRUE, recursive = TRUE)
  ref.file = list.files(paste0(dirs, 'datasets/', data_name[id]), pattern = '._ref.rds', full.names = TRUE, recursive = TRUE)
  Bis.file = list.files(paste0(dirs, 'Results/Bis/', data_name[id]), pattern = '.mat', full.names = TRUE, recursive = TRUE)
  ALRA.file = list.files(paste0(dirs, 'Results/ALRA/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  DCA.file = list.files(paste0(dirs, 'Results/DCA/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  DeepImpute.file = list.files(paste0(dirs, 'Results/DeepImpute/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  MAGIC.file = list.files(paste0(dirs, 'Results/MAGIC/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  SAVER.file = list.files(paste0(dirs, 'Results/SAVER/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  scImpute.file = list.files(paste0(dirs, 'Results/scImpute/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  SCRABBLE.file = list.files(paste0(dirs, 'Results/SCRABBLE/', data_name[id]), pattern = '.mat', full.names = TRUE, recursive = TRUE)
  scScope.file = list.files(paste0(dirs, 'Results/scScope/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  scVI.file = list.files(paste0(dirs, 'Results/scVI/', data_name[id]), pattern = '.rds', full.names = TRUE, recursive = TRUE)
  
  # load data
  dat <- vector("list", 12)
  names(dat) <- c("Dropout", "Ref", "Bis", "ALRA", "DCA", "DeepImpute", "MAGIC", "SAVER", "scImpute", "SCRABBLE", "scScope", "scVI")
  
  for(j in 1:1){
    
    dat[[1]][[j]] = edgeR::cpm(readRDS(samp.file[1]))
    dat[[2]][[j]] = edgeR::cpm(readRDS(ref.file[1]))
    dat[[3]][[j]] = readMat(Bis.file[j])[[1]]
    # dat[[3]][[j]] = as.matrix(read.table(Bis.file[j], sep = ",", header = FALSE))
    dat[[4]][[j]] = readRDS(ALRA.file[j])
    dat[[5]][[j]] = as.matrix(readRDS(DCA.file[j]))
    dat[[6]][[j]] = t(readRDS(DeepImpute.file[j]))
    dat[[7]][[j]] = t(readRDS(MAGIC.file[j]))
    dat[[8]][[j]] = readRDS(SAVER.file[j])
    dat[[9]][[j]] = readRDS(scImpute.file[j])
    dat[[10]][[j]] = readMat(SCRABBLE.file[j])[[1]]
    dat[[11]][[j]] = t(readRDS(scScope.file[j]))
    dat[[12]][[j]] = t(readRDS(scVI.file[j]))
  }
  
  # delete 0
  # for(i in c(1:10)){
  #   index = rowMeans(dat[[1]][[i]]) > 0
  #   for(j in c(1:11)){
  #     dat[[j]][[i]] = dat[[j]][[i]][index, ]
  #   }
  # }
  
  cor.dat <- vector("list", 2)
  names(cor.dat) <- c("gene", "cell")
  
  for (i in c(1:2)) {
    cor.dat[[i]] <- vector("list", 11)
    names(cor.dat[[i]]) <- c("Dropout", "Bis", "ALRA", "DCA", "DeepImpute", "MAGIC", "SAVER", "scImpute", "SCRABBLE", "scScope", "scVI")
    for (j in c(1:10)) {
      cor.dat[[i]][[j]] <- vector("list", 10)
      # names(cor.dat[[i]][[j]]) <- c("drop", "Bis", "ALRA", "DCA", "DeepImpute", "MAGIC", "SAVER", "SCRABBLE", "scScope", "scVI")
    }
  }
  
  for (i in 1:1) {
    
    for (j in 1:11) {
      ind <- c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
      
      cor.dat[[1]][[j]][[i]] <- get.cor.gene(dat[[2]][[i]], dat[[ind[j]]][[i]])
      cor.dat[[2]][[j]][[i]] <- get.cor.cell(dat[[2]][[i]], dat[[ind[j]]][[i]])
      
    }
  }
  
  cor.gene = data.frame(matrix(unlist(cor.dat[[1]]), byrow = T, nrow = 11))
  cor.cell = data.frame(matrix(unlist(cor.dat[[2]]), byrow = T, nrow = 11))
  
  
  
  p.gene = plot_comparison(cor.gene, "Pearson gene-wise correlation", color, dropout_rate = id)
  p.cell = plot_comparison(cor.cell, "Pearson cell-wise correlation", color, dropout_rate = id)
  
  save(file = paste0(dirs, 'Results/Bis/', data_name[id], '/', 'Cor.RData'), cor.gene, cor.cell)
}
