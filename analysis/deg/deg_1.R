library(clusterProfiler)
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
library(org.Hs.eg.db)

my_list=readRDS('/Users/suyc/Project/Bis/results/deg/deg_gene.rds')

bis_gene = setdiff(my_list$our, unlist(my_list[-4]))

write.table(bis_gene, file = '/Users/suyc/Project/Bis/results/deg/bis_1184_gene.csv', sep = ',', 
            row.names = FALSE, col.names = FALSE)

# perform go and kegg
#指定富集分析的物种库
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html


#gene ID转换
gene <- bitr(bis_gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
write.table(gene$ENSEMBL, file = '/Users/suyc/Project/Bis/results/deg/bis_1184_gene_ENSEMBL.csv', sep = ',', row.names = FALSE, col.names = FALSE)

GO<-enrichGO(gene$ENTREZID,#GO富集分析
             OrgDb = GO_database,
             keyType = "ENTREZID",#设定读取的gene ID类型
             ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
             pAdjustMethod = "BH",
             pvalueCutoff  = 1,
             qvalueCutoff  = 1,
             minGSSize = 10, maxGSSize = 5000,
             readable = T)

filer_GO = GO[GO$pvalue < 0.05, asis=T]

bp_GO<-enrichGO(gene$ENTREZID,#GO富集分析
             OrgDb = GO_database,
             keyType = "ENTREZID",#设定读取的gene ID类型
             ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
             pAdjustMethod = "BH",
             pvalueCutoff  = 1,
             qvalueCutoff  = 1,
             minGSSize = 10, maxGSSize = 5000,
             readable = T)

filer_bp_GO = bp_GO[bp_GO$pvalue < 0.05, asis=T]
pdf('/Users/suyc/Project/Bis/plot/deg/revision/deg_mf_graph.pdf', height = 5,width = 5)
plotGOgraph(filer_bp_GO,firstSigNodes=20)
dev.off()


write.table(filer_GO@result, file = "/Users/suyc/Project/Bis/results/deg/revision/go.csv", #将所有GO富集到的基因集所对应的类型写入本地文件从而得到BP/CC/MF
            append = FALSE, quote = TRUE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

R.utils::setOption("clusterProfiler.download.method",'auto')
KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pAdjustMethod = "BH",
                 minGSSize = 10, maxGSSize = 5000,
                 pvalueCutoff = 1,
                 qvalueCutoff = 1)
KEGG <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

KEGG = KEGG[KEGG$pvalue < 0.05, asis=T]

write.table(KEGG@result, file = "/Users/suyc/Project/Bis/results/deg/revision/kegg.csv", #将所有kegg富集到的基因集所对应的类型写入本地文件
            append = FALSE, quote = TRUE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

