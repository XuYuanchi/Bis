source('D:/project/dab/scripts/deg/c_lib.R')
library(patchwork)
library(scales)
library(DESeq2)
library(ggpubr)
library(UpSetR)
color = readRDS("D:/project/Demo_Bulk/color.rds")
a=c("#FED439E5", "#709AE1E5", "#8A9197E5", "#D2AF81E5", "#FD7446E5", "#D5E4A2E5", "#197EC0E5", "#F05C3BE5",
"#46732EE5", "#71D0F5E5", "#370335E5", "#075149E5", "#C80813E5", "#91331FE5", "#1A9993E5", "#FD8CC1E5")
pdf('D:/project/dab/plot/deg/ven.pdf', width = 12, height = 7)
upset(fromList(list(Bulk=blk.deg.genes , 
                    Raw=raw.deg.genes , 
                    Bis=our.deg.genes, 
                    SCRABBLE=SCRABBLE.deg.genes)), nsets = 4, point.size = 2.2, line.size = 0.7, keep.order = T, decreasing = FALSE,
      order.by = c("degree"), set_size.show = T, sets.x.label = "Number of DEGs", mb.ratio = c(0.55, 0.45),
      sets = c("Bulk",'Raw', "Bis", "SCRABBLE"),
      sets.bar.color = "#56B4E9",
      scale.intersections = "identity", scale.sets = "identity",
      intersections = list(list("Bulk"),
                           list("Raw"),
                           list("Bis"),
                           list("SCRABBLE"),
                           list("Bulk",'Raw', 'Bis', "SCRABBLE"),
                           list("Bulk",'Raw', 'Bis'),
                           list("Raw",'Bis', "SCRABBLE"),
                           list("Bulk","Raw"),
                           list("Bulk","Bis"),
                           list("Bulk", "SCRABBLE")),
      queries = list(
        list(query = intersects,params = list("Bulk","Raw"),color = "#709AE1E5",active = T),
        list(query = intersects,params = list("Bulk","Bis"),color = "#709AE1E5",active = T),
        list(query = intersects,params = list("Bulk","SCRABBLE"),color = "#709AE1E5",active = T),
        list(query = intersects,params = list("Bulk"),color = "#FD7446E5",active = T),
        list(query = intersects,params = list("Raw"),color = "#FD7446E5",active = T),
        list(query = intersects,params = list("Bis"),color = "#FD7446E5",active = T),
        list(query = intersects,params = list("SCRABBLE"),color = "#FD7446E5",active = T),
        list(query = intersects,params = list("Bulk",'Raw', 'Bis'),color = "#46732EE5",active = T),
        list(query = intersects,params = list("Raw",'Bis', "SCRABBLE"),color = "#FED439E5",active = T),
        list(query = intersects,params = list("Bulk",'Raw', 'Bis', "SCRABBLE")
             ,color = "#FD8CC1E5",active = T)
      )
)
dev.off()


library(ggplot2)

# 创建数据框
data <- data.frame(
  Method = c("Dropout", "Bis", "SCRABBLE"),
  F.score = c(0.48199446, 0.589446589, 0.536878407),
  AUC = c(0.868069845, 0.882928619, 0.875239751),
  ACC = c(0.764989265, 0.799549667, 0.782007645)
)

data$Method=factor(data$Method, levels = c("Dropout", "Bis", "SCRABBLE"))
# 将数据从宽格式转换为长格式
data_long <- tidyr::gather(data, Metric, Value, -Method)

# 绘制水平的分组柱状图
p=ggplot(data_long, aes(x = Value, y = Metric, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparison of F.score, AUC, and ACC",
       x = "Score", y = "Metric") +
  scale_fill_manual(values = c("#F05C3BE5", "#709AE1E5", "#D5E4A2E5")) +
  coord_cartesian(xlim=c(0.45,0.9))+
  theme_minimal()
ggsave(plot = p, filename = 'D:/project/dab/plot/deg/metric.pdf', height = 5, width = 7)

library(ggpubr)
p=ggbarplot(data_long, "Metric", "Value", fill = "Method", color = "Method",
            palette = c("#FD7446E5", "#709AE1E5", "#D5E4A2E5"), position = position_dodge(0.8))
p=ggpar(p, orientation="horiz")
# p=ggpar(p, ylim = c(0.45,0.9))
ggsave(plot = p, filename = 'D:/project/dab/plot/deg/metric.pdf', height = 5, width = 7)

library(tidyverse)
## plot the correlation scatter plot of log fold changes
plot_correlation <- function(object1, object2, xlab, ylab){
  ovlp <- intersect(rownames(object1[object1$padj <= 0.05,]),rownames(object2[object2$padj <= 0.05,]))
  x = object1[ovlp,]$log2FoldChange
  y = object2[ovlp,]$log2FoldChange
  cor <- cor.test(x,y)
  r <- round(cor$estimate,2) 
  p <- signif(cor$p.value,3)
  up <- sum(x>0)
  down <- sum(x<0)
  data <- data.frame(x=x, y=y)
  p <- data %>% mutate(Color = ifelse(x*y>0, ifelse(x+y > 0, "#C80813E5", "#709AE1E5"), "gray")) %>%
    ggplot(aes(x=x, y=y, color = Color)) +
    geom_point() +
    stat_smooth(method = 'lm', colour = '#FD7446E5') +
    annotate("text", x = -Inf, y = Inf, label = paste0("r = ", r), hjust = -0.2, 
             vjust = 2, size = 5) +
    annotate("text", x = -Inf, y = Inf, label = paste0("up = ", up), hjust = -0.2, 
             vjust = 4, size = 5) +
    annotate("text", x = -Inf, y = Inf, label = paste0("down = ", down), hjust = -0.2, 
             vjust = 6, size = 5) +
    scale_color_identity() +
    xlab(NULL) +
    theme_bw() +
    labs(x = xlab, y = ylab) +
    theme(
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 12),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  return(p)
}



p1=plot_correlation(blk.result[blk.top1000.genes,], raw.result[rownames(raw.result) %in% blk.top1000.genes,], "Bulk", "Raw")
ggsave(plot=p, filename = 'D:/project/dab/plot/deg/cor_raw.pdf', height = 5, width = 5)
p2=plot_correlation(blk.result[blk.top1000.genes,], scimpute.result[blk.top1000.genes,], "Bulk", "Bis")
ggsave(plot=p, filename = 'D:/project/dab/plot/deg/cor_bis.pdf', height = 5, width = 5)
p3=plot_correlation(blk.result[blk.top1000.genes,], SCRABBLE.result[rownames(SCRABBLE.result) %in% blk.top1000.genes,], "Bulk", "SCRABBLE")
ggsave(plot=p, filename = 'D:/project/dab/plot/deg/cor_scrabble.pdf', height = 5, width = 5)
ggsave(plot = p1+p2+p3+plot_layout(design = "ABC"), filename = 'D:/project/dab/plot/deg/cor.pdf', height = 5, width = 15)

plot.gene <- function (dds, gene, intgroup = "condition", normalized = TRUE, 
                       transform = F, main, xlab = "group", returnData = TRUE, 
                       replaced = FALSE, pc, mar = c(2,4,4,0.3), ...) 
{
  set.seed(200)
  stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) & 
                                                         (gene >= 1 & gene <= nrow(dds)))))
  if (!all(intgroup %in% names(dds@colData))) 
    stop("all variables in 'intgroup' must be columns of colData")
  if (!returnData) {
    if (!all(sapply(intgroup, function(v) is(dds@colData[[v]], 
                                             "factor")))) {
      stop("all variables in 'intgroup' should be factors, or choose returnData=TRUE and plot manually")
    }
  }
  if (missing(pc)) {
    pc <- if (transform) 
      0.5
    else 0
  }
  if (is.null(dds@colData@listData[["sizeFactor"]]) & is.null(normalizationFactors(dds))) {
    dds <- estimateSizeFactors(dds)
  }
  cnts <- counts(dds, normalized = normalized, replaced = replaced)[gene, ]
  group <- if (length(intgroup) == 1) {
    dds@colData[[intgroup]]
  }
  else if (length(intgroup) == 2) {
    lvls <- as.vector(t(outer(levels(dds@colData[[intgroup[1]]]), 
                              levels(dds@colData[[intgroup[2]]]), function(x, 
                                                                           y) paste(x, y, sep = ":"))))
    droplevels(factor(apply(as.data.frame(dds@colData[, 
                                                      intgroup, drop = FALSE]), 1, paste, collapse = ":"), 
                      levels = lvls))
  }
  else {
    factor(apply(as.data.frame(dds@colData[, intgroup, 
                                           drop = FALSE]), 1, paste, collapse = ":"))
  }
  data <- data.frame(count = cnts + pc, group = as.integer(group))
  logxy <- if (transform) 
    "y"
  else ""
  if (missing(main)) {
    main <- if (is.numeric(gene)) {
      rownames(dds)[gene]
    }
    else {
      gene
    }
  }
  ylab <- ifelse(normalized, "normalized count", "count")
  if (returnData) 
    return(data.frame(count = data$count, dds@colData[intgroup]))
  # par(mar = mar)
  # boxplot(count ~ group, data = data, alpha = 0.4, outline = F,
  #         xlim = c(0.5, max(data$group) + 0.5), log = logxy, xaxt = "n", 
  #         xlab = xlab, ylab = ylab, main = main, ...)
  # rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
  #        "white")
  # boxplot(count ~ group, data = data, alpha = 0.4, outline = F,
  #         xlim = c(0.5, max(data$group) + 0.5), log = logxy, xaxt = "n", 
  #         xlab = xlab, ylab = ylab, main = main, add = T, ...)
  # points( jitter(data$group) , data$count, pch = 16, col = alpha(c("#4682B466", "#FFBB0066")[data$group], 0.4) )
  # axis(1, at = seq_along(levels(group)), levels(group))
  
  # 创建小提琴图
  # ggplot(data, aes(x = group, y = count)) +
  #   geom_violin(alpha = 0.4, fill = "#4682B466") +
  #   geom_jitter(width = 0.1, alpha = 0.4, color = "#FFBB0066") +
  #   xlab(xlab) +
  #   ylab(ylab) +
  #   ggtitle(main) +
  #   theme_minimal()
  
}

gene = c('NANOG', 'SOX2', 'DNMT3B', 'POU5F1', 'ZFP42', 'GATA6', 'CER1', 'EOMES', 'LEFTY1', 'CXCR4')
for(i in c(1:10)){
  p=plot.gene(SCRABBLE.deg, gene[i], intgroup = 'Cell', xlab = "n", main = "SCRABBLE",  mar = c(0,2,2,0.3))
  
  fig=ggboxplot(p, x="Cell", y="count", add = c("jitter"), fill = "Cell", title = 'SCRABBLE', 
                palette = c("#1A9993E5", "#E7B800"), add.params = list(fill = "white", size=0.1)) + border("black")
  fig = ggpar(fig, legend = "none")
  fig = annotate_figure(fig, right = gene[i])
  
  ggsave(plot=fig, filename = paste0('D:/project/dab/plot/deg/scrabble_', gene[i],'.pdf'), height = 3.5, width = 4)
}
## plot 
{
  pdf('D:/project/dab/plot/deg/gene_NANOG.pdf', width = 3.6, height = 2.5)
  par(mfrow=c(1,3),cex = 0.6, bg = "white")
  plot.gene(blk.deg, "NANOG", intgroup = 'Cell', xlab = "n", main = "Bulk",  mar = c(0,2,2,0.3))
  plot.gene(raw.deg, "NANOG", intgroup = 'Cell', xlab = "n", main = "Raw",  mar = c(0,2,2,0.3))
  plot.gene(our.deg, "NANOG", intgroup = 'Cell', xlab = "n", main = "Bis",  mar = c(0,2,2,0.3))
  plot.gene(SCRABBLE.deg, "NANOG", intgroup = 'Cell', xlab = "n", main = "SCRABBLE",  mar = c(0,2,2,0.3))
  dev.off()
}


bis_raw_diff = setdiff(our.deg.genes,raw.deg.genes)
bis_scrabble_diff = setdiff(our.deg.genes,SCRABBLE.deg.genes)
bis_diff = intersect(bis_raw_diff, bis_scrabble_diff)
bis_bulk_inter = intersect(bis_diff, blk.deg.genes)

write.table(bis_bulk_inter, file = 'D:/project/dab/plot/deg/deg.csv', sep = ',')
