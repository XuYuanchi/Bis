library(monocle3)
library(gridExtra)
library(ggplot2)
library(scales)
source("analysis_lib.R")

cell_meta <- read.table("data/Timecourse.label.txt", header = F)

## raw
d <- read.table(file = "data/sc/sc_raw.csv", sep = ",", header = T)

gene_short_name <- d[, 1]

rownames(d) <- gene_short_name

d <- data.matrix(d[, -1])

rownames(cell_meta) <- colnames(d)
colnames(cell_meta) <- "Time"

gene.anno <- data.frame(gene_short_name, row.names = gene_short_name)

data_raw <- generate_ICT(d, cell_meta, gene.anno, "raw")

## SCRABBLE
d <- read.table(file = "Results/SCRABBLE_impute.csv", sep = ",")
rownames(d) <- gene_short_name
colnames(d) <- rownames(cell_meta)
d <- data.matrix(d)
data_SCRABBLE <- generate_ICT(d, cell_meta, gene.anno, "SCRABBLE")
## DCA
d <- read.table(file = "Results/DCA_impute.tsv", sep = "\t", header = T)
rownames(d) <- d[, 1]
d <- data.matrix(d[, -1])
data_DCA <- generate_ICT(d, cell_meta, gene.anno, "DCA")
## SAVER
d <- read.table(file = "Results/SAVER_impute.csv", sep = ",")
rownames(d) <- gene_short_name
colnames(d) <- rownames(cell_meta)
d <- data.matrix(d)
data_SAVER <- generate_ICT(d, cell_meta, gene.anno, "SAVER")
## scImpute
d <- read.table(file = "Results/scImpute_impute.csv", sep = ",", header = T)
rownames(d) <- d[, 1]
d <- data.matrix(d[, -1])
data_scImpute <- generate_ICT(d, cell_meta, gene.anno, "scImpute")
## DeepImpute
d <- read.table(file = "Results/deepImpute.csv", sep = ",", header = T)
rownames(d) <- d[, 1]
d <- data.matrix(d[, -1])
data_DeepImpute <- generate_ICT(d, cell_meta, gene.anno, "DeepImpute")
## MAGIC
d <- read.table(file = "Results/MAGIC_impute.csv", sep = ",")
rownames(d) <- gene_short_name
colnames(d) <- rownames(cell_meta)
d <- data.matrix(d)
data_MAGIC <- generate_ICT(d, cell_meta, gene.anno, "MAGIC")
## Bis
d <- read.table(file = "Results/Bis_impute.csv", sep = ",")
rownames(d) <- gene_short_name
colnames(d) <- rownames(cell_meta)
d <- data.matrix(d)
data_Bis <- generate_ICT(d, cell_meta, gene.anno, "Bis")

save(data_raw, data_SCRABBLE, data_DCA, data_SAVER, data_scImpute, data_DeepImpute, data_MAGIC, data_Bis, file = "Results/Timecourse_Monocle_data.Rda")
