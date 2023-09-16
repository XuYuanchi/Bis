# boxplot
library(ggpubr)

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

p=ggbarplot(data_long, x='Metric', y='Value', palette = c("#FED439E5", "#709AE1E5", "#46732EE5"),
          fill = 'Method', position = position_dodge(0.9), color='black')

ggpar(p, ylim = c(0.45,0.9))

ggsave(filename="/Users/suyc/Project/dab/plot/deg/metric_new.pdf", height = 4, width =4, bg = "white")
