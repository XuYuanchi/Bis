# plot tracjectory inference figure
library(patchwork)
library(ggplot2)
library(ggpubr)
library(tidyr)
color = c("#F05C3BE5", "#FED439E5", "#46732EE5", "#D5E4A2E5",  "#709AE1E5", "#FD8CC1E5")
color1 = c("#FED439E5", "#709AE1E5", "#8A9197E5", "#D2AF81E5", "#FD7446E5", "#D5E4A2E5", "#197EC0E5", "#F05C3BE5",
           "#46732EE5", "#71D0F5E5", "#370335E5")
# load data
re =  readRDS('results/ti/re.rds')

p = re[[2]][[1]]
for(i in c(3:11)){
  p = p + re[[i]][[1]] 
}
fig_layout <- "
ABCDE
FGHIJ
"

p = p + plot_layout(design = fig_layout, guides = 'collect') & 
  theme(legend.position = 'bottom', panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid")) & 
  scale_color_manual(values = color) 

ggsave(plot=p, filename = 'plot/ti/trac.pdf', height = 5, width = 10)
ggsave(plot = re[[1]][[1]]& 
         theme(legend.position = 'bottom', panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid")) & 
         scale_color_manual(values = color), filename = 'plot/ti/trac_raw.pdf', height = 5, width = 4)

# plot metric 
df = data.frame(Methods = c('Dropout', 'Bis', 'ALRA', 'DCA', 'DeepImpute', 'MAGIC', 'SAVER', 
                              'scImpute', 'SCRABBLE', 'scScope', 'scVI'))

for(i in c(1:11)){
  df[i,2] = re[[i]][[2]]
  df[i,3] = re[[i]][[3]]
}
colnames(df) = c('Methods', 'cor', 'POS')
# Pseudo-temporal Ordering(POS)
# Kendall’s rank correlation scores

# ggplot(df, aes(Methods, cor, fill=Methods)) +
#   geom_bar(stat = "identity") + 
#   scale_fill_manual(values = color1) + 
#   labs(title = "a", x = "methods", y = "value")

p1 = ggbarplot(df, x = "Methods", y = "cor", fill = "Methods", palette = color1, 
              width = 0.8, ylab = "Kendall's rank correlation", ggtheme = theme_bw()) + 
  rotate_x_text(45)

p2 = ggbarplot(df, x = "Methods", y = "POS", fill = "Methods", palette = color1, 
          width = 0.8, ylab = "Pseudo-temporal Ordering(POS)", ggtheme = theme_bw()) + 
  rotate_x_text(45)
p = p1 + p2 + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
ggsave(plot = p, filename = 'plot/ti/trac_metric.pdf', height = 6, width = 9)


