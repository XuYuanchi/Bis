library(monocle3)
library(TSCAN)
library(ggplot2)
library(R.matlab)
library(scales)
library(rhdf5)

source('scripts/ti/utils.R')

re = list()
# dropout
dat = read.table(file = 'data/ti/gse.csv', sep = ',', header = T)
rownames(dat) = dat[,1]
dat = dat[,-1]

gene_short_name = row.names(dat)
cell_meta = read.table(file = 'data/ti/label.txt', header = F)
rownames(cell_meta) = colnames(dat)
colnames(cell_meta) = "Time"
cell_label = factor(cell_meta$Time, levels = unique(cell_meta$Time))

gene.anno <- data.frame(gene_short_name, row.names = gene_short_name)

re[[1]] = process(as.matrix(dat), cell_meta, gene.anno, 'Dropout', cell_label, levels(cell_label)[1])

# Bis
dat = readMat('results/ti/gse_1.mat')[[1]]
re[[2]] = process(dat, cell_meta, gene.anno, 'Bis', cell_label, levels(cell_label)[1])
re[[2]][[1]]
re[[2]][[2]]
re[[2]][[3]]

# ALRA
dat = readRDS('results/ti/alra_gse.rds')
re[[3]] = process(dat, cell_meta, gene.anno, 'ALRA', cell_label, levels(cell_label)[1])

# dca
dat = read.table('results/ti/dca_gse.tsv', sep = '\t', header = T)
dat = readMat('results/ti/dca_gse.mat')[[1]]
# rownames(dat) = dat[,1]
# dat = dat[,-1]
# gene_short_name1 = row.names(dat)
# gene.anno1 <- data.frame(gene_short_name1, row.names = gene_short_name1)
re[[4]] = process(dat, cell_meta, gene.anno, 'DCA', cell_label, levels(cell_label)[1])

# DeepImpute
#dat = readMat('results/ti/deepimpute_gse.mat')[[1]]
h5file=H5Fopen('results/ti/deepimpute_gse.h5')
dat = h5file$data 
re[[5]] = process(t(dat), cell_meta, gene.anno, 'DeepImpute', cell_label, levels(cell_label)[1])

# MAGIC
dat = readMat('results/ti/magic_gse.mat')[[1]]
re[[6]] = process(dat, cell_meta, gene.anno, 'MAGIC', cell_label, levels(cell_label)[1])

# SAVER
dat = readRDS('results/ti/saver_gse.rds')
re[[7]] = process(dat, cell_meta, gene.anno, 'SAVER', cell_label, levels(cell_label)[1])

# scImpute
dat = readRDS('results/ti/scimpute_gse.rds')
re[[8]] = process(dat, cell_meta, gene.anno, 'scImpute', cell_label, levels(cell_label)[1])

# SCRABBLE
dat = readMat('results/ti/scrabble_gse.mat')[[1]]
dat = read.table('results/ti/scrabble_gse.csv', sep = ',', header = F)
colnames(dat) = rownames(cell_meta)
re[[9]] = process(as.matrix(dat), cell_meta, gene.anno, 'SCRABBLE', cell_label, levels(cell_label)[1])

# scScope
dat = readMat('results/ti/scscope_gse.mat')[[1]]
re[[10]] = process(dat, cell_meta, gene.anno, 'scScope', cell_label, levels(cell_label)[1])

# scVI
dat = readMat('results/ti/scvi_gse.mat')[[1]]
re[[11]] = process(dat, cell_meta, gene.anno, 'scVI', cell_label, levels(cell_label)[1])

saveRDS(re, file='results/ti/re.rds')
