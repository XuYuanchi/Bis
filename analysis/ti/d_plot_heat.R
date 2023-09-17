# Differential gene expression usig monocle3
library(monocle3)
library(R.matlab)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(grid)
source('scripts/ti/utils.R')
re = readRDS('results/ti/re.rds')
mycds2 = re[[2]][[4]]
mycds2_res <- graph_test(mycds2, neighbor_graph="principal_graph", cores=4)#拟时差异基因
genes1 <- row.names(subset(mycds2_res, q_value< 0.01 & morans_I > 0.2))#我们选择显著的基因进行

mycds2 = re[[1]][[4]]
mycds2_res <- graph_test(mycds2, neighbor_graph="principal_graph", cores=4)#拟时差异基因
genes2 <- row.names(subset(mycds2_res, q_value< 0.01 & morans_I > 0.2))#我们选择显著的基因进行

genes = setdiff(genes1, genes2)

raw_dat = read.table(file = 'data/ti/gse.csv', sep = ',', header = T)
rownames(raw_dat) = raw_dat[,1]
raw_dat = raw_dat[,-1]

dat = readMat('results/ti/gse_1.mat')[[1]]
rownames(dat) = rownames(raw_dat)
colnames(dat) = colnames(raw_dat)
mycds2 = re[[2]][[4]]

plot_matrix <- dat[match(genes,#挑选我们前面确定需要展示的基因
                                   rownames(rowData(mycds2))),
                             order(pseudotime(mycds2))]#按照拟时对基因排序

#数据拟合平滑处理smooth.spline。计算Z-score
plot_matrix <- t(apply(plot_matrix,1,function(x){smooth.spline(x,df=3)$y}))
plot_matrix <- t(apply(plot_matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(plot_matrix) <- genes;
dim(plot_matrix)


#得到矩阵就可以做热图了，很简单
##排序设置
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}


cutColumn_Means <- function(data_exp,#需要分割数据
                            cut#需要分割的列数
){
  plot_matrix_combin <- list()
  nums <- ncol(data_exp)/cut
  if (nums-round(nums, 0)==0){
    
    for (i in 1:length(seq(1, ncol(data_exp), cut))){
      num <- seq(1, ncol(data_exp), cut)
      A <- as.data.frame(rowMeans(data_exp[,num[i]:(cut+num[i]-1)]))[,1]
      plot_matrix_combin[[i]] <- A
      
    }
    plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
    rownames(plot_matrix_combin) <- rownames(data_exp)
    colnames(plot_matrix_combin) <- seq(1,ncol(plot_matrix_combin),1)
    return(plot_matrix_combin)
    
  }else{
    
    for (i in 1:length(seq(1, ncol(data_exp)-cut, cut))){
      num <- seq(1, ncol(data_exp)-cut, cut)
      A <- as.data.frame(rowMeans(data_exp[,num[i]:(cut+num[i]-1)]))[,1]
      plot_matrix_combin[[i]] <- A
    }
    
    plot_matrix_combin[[length(seq(1, ncol(data_exp)-cut, cut))+1]] <- as.data.frame(rowMeans(data_exp[,(max(num)+cut):ncol(data_exp)]))                       
    plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
    rownames(plot_matrix_combin) <- rownames(data_exp)
    colnames(plot_matrix_combin) <- seq(1,ncol(plot_matrix_combin),1)
    return(plot_matrix_combin)
  }
  
  
}

plot_test <- cutColumn_Means(plot_matrix,cut = 50)

# 删掉标准差为0的行
plot_test = plot_test[apply(plot_test, 1, function(x) sd(x)!=0),] 

# 删掉标准差为0的列
plot_test = plot_test[,apply(plot_test, 2, function(x) sd(x)!=0)]

#plot_test[is.na(plot_test)] = 0
plot_test = na.omit(plot_test)

p <- pheatmap::pheatmap(plot_test, 
                         useRaster = T,
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=F, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         cutree_rows=6,
                         filename=NA,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#709AE1E5','white',"#F05C3BE5"))(100),
                         clustering_callback = callback)

#行注释
annotation_row <- data.frame(Cluster=factor(cutree(p$tree_row, 6)))
row.names(annotation_row) <- rownames(plot_test)

#rowcolor <- c("#85B22E","#E29827","#922927",'#57C3F3') 
rowcolor = c("#6BB952", "#EC748B", "#C4A751",
             "#36ACA2", "#6EB1DE", "#B67FB3")
rowcolor = c("#f7c0c5", "#ffe5be", "#c5e9c7",
             "#bee8f0", "#bfcfe4", "#e2c5de")
rowcolor = c("#a7cee2", "#b4d88a", "#34a048",
             "#f69999", "#fdbf6e", "#f57f20")
names(rowcolor) <- c("1","2","3","4", "5", "6") #类型颜色
#注释颜色设置
ann_colors <- list(Cluster=rowcolor) #颜色设置


p <- pheatmap::pheatmap(plot_test, 
                         cluster_cols=FALSE, 
                         cluster_rows=T, 
                         show_rownames=F, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         filename=NA,
                         cutree_rows=6,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#2078b4','white',"#e21f26"))(100),
                         annotation_colors=ann_colors,
                         annotation_row = annotation_row,
                         clustering_callback = callback,
                         annotation_names_col = F,
                         annotation_names_row = F,
                         main="Pseudotime")
p


#展示需要的基因
gene <- c('POU5F1', 'NANOG', 'SOX2', 'DNMT3B', 'NODAL', 'EOMES', 'ID1', 'CDX1', 'T', 'MSX2', 'CER1', 'GATA4', 'DKK4',
          'MYCT1', 'POU2AF1', 'PRDM1', 'CXCR4', 'SOX17', 'ACADM', 'MMP9', 'PAF1', 'SETD2', 'SOX17', 
          'ISL1', 'LEF1', 'WNT5A', 'PSEN1', 'RWDD1', 'SRC','TAF1','TSG101',
          'CREBZF', 'EIF2A', 'ERLIN2', 'FBXW7', 'HSPA5', 'INSIG1',  'YAP1', 'TCF7L2', 'TMEM170B', 'TMEM9', 'TNKS', 
          'TSC2', 'UBAC2', 'UBR5', 'VPS35')
gene = rownames(plot_test)[1:50]
gene = c('FOXRED2','GCH1','GEN1','GIGYF2','GPX8','GSTM1','HILPDA','HSPA13','HSPA1A','HTRA2','HUS1','IGFBP6','INPP5F','INSIG1','KAT5','KAT6A','KDM4D','KLF4','KLHDC10',
         'ITPR1','LIMD1','MECP2','NDRG1','NOL3','PDK1','PIK3CB','PINK1','PLOD1','POLB','RORA','SLC1A1','SMAD3','SMAD4','VEGFB',
         'ABL1','ACTL6A','AIDA','AQP11','ATF4','ATXN7L3','AXIN1','AXIN2','BABAM1','BAD','BAG5','BAG6','BAK1','BAX','BCL2L11','BFAR','BMP7','BOK','CCDC117','CEBPG','CHEK2',
         'CUL4A','DDRGK1','DDX11','DHX9','DMAP1','DNAJB1', 'YAP1',
         'ADNP2','ALKBH3','ASCC2','ASCC3','ATG5','ATMIN','ATR','AUP1')
#这里我们使用add.flag.R函数添加，需要注意的是添加flag要成功，不要添加cutree_rows参数
p <- pheatmap::pheatmap(plot_test, 
                         cluster_cols=FALSE, 
                         cluster_rows=T,  
                         show_rownames=T, 
                         show_colnames=F, 
                         clustering_method = "ward.D2",
                         filename=NA,
                         border_color = NA,
                         fontsize_row = 8,
                         color=colorRampPalette(c('#2078b4','white',"#e21f26"))(100),
                         annotation_colors=ann_colors,
                         annotation_row = annotation_row,
                         clustering_callback = callback,
                         annotation_names_col = F,
                         annotation_names_row = F,
                         main="Pseudotime")

p = add.flag(p,kept.labels = gene,repel.degree = 0.2)
ggsave(plot = p, filename = 'plot/ti/heat.pdf', height = 10, width = 10)

