# cal acc..
library(ROCR)
precision.recall <- function(c1, c2, plot = F, label = "ROC.curve") {
  # ref.label <- as.vector(outer(c1, c1, "=="))
  # pred.label <- as.vector(outer(c2, c2, "=="))
  # pred <- prediction(as.numeric(pred.label), as.numeric(ref.label) )
  pred <- prediction(as.numeric(c1), as.numeric(c2))

  # Accuracy
  acc.tmp <- performance(pred, "acc")
  acc <- as.numeric(acc.tmp@y.values[[1]][2])

  # ROC curve
  ROC.perf <- performance(pred, "tpr", "fpr")
  # ROC area under the curve
  auc.tmp <- performance(pred, "auc")
  auc <- as.numeric(auc.tmp@y.values)

  if (plot) {
    # pdf(paste(label,".pdf", sep = ""),2.5, 3.0)
    plot(ROC.perf)
    text(0.6, 0.05, paste("AUC=", round(auc, digits = 4)))
    abline(a = 0, b = 1)
    # dev.off()
  }

  ## precision
  prec.tmp <- performance(pred, "ppv")
  prec <- as.numeric(prec.tmp@y.values[[1]][2])
  ## F1-score
  f.tmp <- performance(pred, "f")
  f <- as.numeric(f.tmp@y.values[[1]][2])
  ##
  return(list(F.score = f, AUC = auc, ACC = acc))
}

load("/data/DEG.DESeq2.results.rd")
load("D:/project/dab/results/deg/DEG_our.rd")

lfc <- 1.5

blk.deg.genes <- rownames(blk.result[blk.result$padj <= 0.05 & abs(blk.result$log2FoldChange) > lfc & blk.result$baseMean >= 10, ])
our.deg.genes <- rownames(our.result[our.result$padj <= 0.05 & abs(our.result$log2FoldChange) > lfc & our.result$baseMean >= 10, ])
SCRABBLE.deg.genes <- rownames(SCRABBLE.result[SCRABBLE.result$padj <= 0.05 & abs(SCRABBLE.result$log2FoldChange) > lfc & SCRABBLE.result$baseMean >= 10, ])

all.degs <- unique(c(
  blk.deg.genes,
  raw.deg.genes,
  our.deg.genes,
  SCRABBLE.deg.genes
))

all.gene <- rownames(blk.cts)

blk.deg.pred <- all.gene %in% all.degs
raw.deg.pred <- all.gene %in% raw.deg.genes
our.deg.pred <- all.gene %in% our.deg.genes
SCRABBLE.deg.pred <- all.gene %in% SCRABBLE.deg.genes

DEG.performance <- rbind(
  unlist(precision.recall(blk.deg.pred, raw.deg.pred)),
  unlist(precision.recall(blk.deg.pred, our.deg.pred)),
  unlist(precision.recall(blk.deg.pred, SCRABBLE.deg.pred))
)
DEG.performance

rownames(DEG.performance) <- c("Raw", "Bis", "SCRABBLE")

save(DEG.performance, file = "D:/project/dab/results/DEG.performance.RData")

write.table(DEG.performance, file = "D:/project/dab/results/DEG.performance.csv", sep = ",", col.names = NA)
