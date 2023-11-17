# plot bubble figure
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

##气泡图，多组情况
go_enrich <- read.table('/Users/suyc/Project/Bis/results/deg/revision/filter_go.csv', sep = ',',
                        header = T)
go_enrich = go_enrich[,-1]


fenshu2xiaoshu<-function(ratio){
  sapply(ratio,function(x) as.numeric(strsplit(x,"/")[[1]][1])/as.numeric(strsplit(x,"/")[[1]][2]))
}
go_enrich$geneRatio=fenshu2xiaoshu(go_enrich$GeneRatio)

#纵坐标是 GO Term，横坐标是各个分组比较
#按各 GO Term 中富集的基因数量（Count）赋值气泡图中点的大小，颜色按 p 值着色
ggplot(go_enrich[go_enrich$ONTOLOGY=='BP', ], aes(geneRatio, Description)) +
  geom_point(aes(color = -log10(pvalue), size = Count)) +
  scale_size(range = c(2, 6)) +
  scale_color_gradientn(colors = c("#57121d", "#d56e5e", "#eaebea", "#5390b5", "#1f294e")) +
  #scale_y_discrete(labels = function(x) str_wrap(x, width = 60) ) +
  theme_bw() +
  ggtitle('Biological Process') +
  labs(x = 'GeneRatio', y = '')
ggsave(filename = '/Users/suyc/Project/Bis/plot/deg/revision/bp.pdf', height = 8, width =8)

#纵坐标是 GO Term，横坐标是各个分组比较
#按各 GO Term 中富集的基因数量（Count）赋值气泡图中点的大小，颜色按 p 值着色
ggplot(go_enrich[go_enrich$ONTOLOGY=='CC', ], aes(geneRatio, Description)) +
  geom_point(aes(color = -log10(pvalue), size = Count)) +
  scale_size(range = c(2, 6)) +
  scale_color_gradientn(colors = c("#57121d", "#d56e5e", "#eaebea", "#5390b5", "#1f294e")) +
  #scale_y_discrete(labels = function(x) str_wrap(x, width = 60) ) +
  theme_bw() +
  ggtitle('Cellular Component') +
  labs(x = 'GeneRatio', y = '')
ggsave(filename = '/Users/suyc/Project/Bis/plot/deg/revision/cc.pdf', height = 8, width =8)

#纵坐标是 GO Term，横坐标是各个分组比较
#按各 GO Term 中富集的基因数量（Count）赋值气泡图中点的大小，颜色按 p 值着色
ggplot(go_enrich[go_enrich$ONTOLOGY=='MF', ], aes(geneRatio, Description)) +
  geom_point(aes(color = -log10(pvalue), size = Count)) +
  scale_size(range = c(2, 6)) +
  scale_color_gradientn(colors = c("#57121d", "#d56e5e", "#eaebea", "#5390b5", "#1f294e")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50) ) +
  theme_bw() +
  ggtitle('Molecular Function') +
  labs(x = 'GeneRatio', y = '')
ggsave(filename = '/Users/suyc/Project/Bis/plot/deg/revision/mf.pdf', height = 8, width =8)


# kegg
#keggrich = read.table('/Users/suyc/Project/Bis/results/deg/revision/filter_go.csv', sep = ',',
                        #header = T)