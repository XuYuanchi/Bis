# analysis using deseq2

library(DESeq2)
library(R.matlab)

# raw 
raw.cts = read.table('D:/project/dab/data/deg/sc_raw.csv', sep = ",", header = T, row.names = 1)
raw.cts <- round(raw.cts)
raw.label <- read.table('D:/project/dab/data/deg/DEG.label.txt')
row.names(raw.label) <- colnames(raw.cts)
colnames(raw.label) <- "Cell"

# load data
our.cts = readMat('D:/project/dab/results/deg_raw.mat')[[1]]

rownames(our.cts) <- rownames(raw.cts)
colnames(our.cts) <- colnames(raw.cts)
our.cts <- round(data.matrix(our.cts))
our.dds <- DESeqDataSetFromMatrix(countData = our.cts,
                                  colData = raw.label,
                                  design = ~ Cell)

our.deg <- DESeq(our.dds)
our.result <- results(our.deg, contrast = c('Cell', 'H1', 'DEC'))
our.result <- our.result[complete.cases(our.result),]
our.deg.list <- our.result[our.result$padj <= 0.05 & abs(our.result$log2FoldChange) >1.5 & our.result$baseMean >= 10,]
our.deg.up.num <- sum(our.deg.list$log2FoldChange>0)
our.deg.down.num <- sum(our.deg.list$log2FoldChange<0)

our.deg.genes <- rownames(our.deg.list)
save(our.deg.genes, our.result, our.deg, file = "D:/project/dab/results/deg/DEG_our.rd")

