# gene go analysis

library(clusterProfiler)
library(org.Hs.eg.db)

###首先提取热图中各个module的基因
module_gene <- as.data.frame(cutree(p$tree_row, k=6))
colnames(module_gene) <- "Module"
module_gene$gene <- rownames(module_gene)
#我们这里展示GO结果
Module_GO=data.frame()

for (i in unique(module_gene$Module)) {
  
  data=filter(module_gene,module_gene$Module==i)
  df=bitr(data$gene, 
          fromType="SYMBOL",
          toType=c("ENTREZID"), 
          OrgDb="org.Hs.eg.db")#Symbol转化为ID
  
  go <- enrichGO(gene= unique(df$ENTREZID),
                 OrgDb= org.Hs.eg.db,
                 keyType= 'ENTREZID',
                 ont= "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff= 0.05,
                 qvalueCutoff= 0.05,
                 readable= TRUE,
                 minGSSize = 10, maxGSSize = 5000)
  go_res=go@result
  
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    Module_GO=rbind(Module_GO,go_res)
  }
}

#筛选显著的Terms
Module_GO <- Module_GO[which(Module_GO$qvalue <= 0.05),]
Module_GO <- Module_GO[,c("ID","Description","qvalue","cluster", "geneID")]
write.csv(Module_GO, file = 'plot/ti/Module_GO_MF.csv')


dat = read.table(file = 'data/batch/uc3.csv', sep = ',')
library(R.matlab)
dataset_name = c("cell_lines", "panc8_rm", "uc3", "crc")
i=4
dat = readMat(paste0('data/batch/', dataset_name[i], '.mat'))[[1]]
saveRDS(t(dat), file = paste0("data/batch/", dataset_name[i], ".rds"))

