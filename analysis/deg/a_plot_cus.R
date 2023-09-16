library(tidyverse) 
library(patchwork) 
library(ComplexHeatmap)

library(RVenn) # Only required if you want Venn diagrams 
library(RColorBrewer) # This is for the colors only, not actually necessary

my_list=readRDS('/Users/suyc/Project/dab/results/deg/deg_gene.rds')
names(my_list) = c('Bulk', 'Raw', 'SCRABBLE', 'Bis')
# my_object <- RVenn::Venn(my_list)
# 
# ggvenn(
#   my_object, slice = 1:3, 
#   thickness = 0.5,
#   alpha = 0.5, 
#   fill = brewer.pal(8, "Set2")
# ) +
#   theme_void() +
#   theme(
#     legend.position = "none"
#   )
# 
# ggsave("../Results/VennDiagram_quick_start.svg", height = 4, width = 4, bg = "white")
# ggsave("../Results/VennDiagram_quick_start.png", height = 4, width = 4, bg = "white")


# ComplexHeatmap for heavy lifting
comb_mat <- make_comb_mat(my_list)
my_names <- set_name(comb_mat)
#my_names <- c('Bulk', 'Raw', 'SCRABBLE', 'Bis')
# Total set size
my_set_sizes <- set_size(comb_mat) %>% 
  as.data.frame() %>% 
  rename(sizes = ".") %>% 
  mutate(Set = row.names(.)) 
color = c("#FD8CC1E5","#FED439E5","#46732EE5", "#709AE1E5")
p1 <- my_set_sizes %>% 
  mutate(Set = reorder(Set, sizes)) %>% 
  ggplot(aes(x = Set, y= sizes)) +
  geom_bar(stat = "identity", aes(fill = Set), alpha = 0.8, width = 0.7, colour='black') +
  geom_text(aes(label = sizes), 
            size = 5, angle = 0, hjust = 0, y = 1) +
  scale_fill_manual(values = color,  # feel free to use some other colors  
                    limits = my_names) + 
  labs(x = NULL,
       y = "Set size",
       fill = NULL) +
  theme_classic() +
  theme(legend.position = "right",
        text = element_text(size= 14),
        axis.ticks.x = element_blank(),
        axis.text = element_blank()
  ) +
  coord_flip()
p1
# Legend
get_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

p2 <- get_legend(p1)
# Overlap sizes
my_overlap_sizes <- comb_size(comb_mat) %>% 
  as.data.frame() %>% 
  rename(overlap_sizes = ".") %>% 
  mutate(category = row.names(.))

p3 <- my_overlap_sizes %>% 
  mutate(category = reorder(category, -overlap_sizes)) %>% 
  ggplot(aes(x = category, y = overlap_sizes)) +
  geom_bar(stat = "identity", fill = "grey80", color = NA, alpha = 0.8, width = 0.7) +
  geom_text(aes(label = overlap_sizes, y = 0), 
            size = 5, hjust = 0, vjust = 0.5) +
  labs(y = "Intersect sizes",
       x = NULL) +
  theme_classic() +
  theme(text = element_text(size= 14, color = "black"),
        axis.text =element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(hjust = 0),
  ) +
  coord_flip()

colors <- colorRampPalette(c("#A9A9A9", "#d3D3D3"))(15)
colors

p3 <- my_overlap_sizes %>% 
  mutate(category = reorder(category, -overlap_sizes)) %>% 
  ggplot(aes(x = category, y = overlap_sizes, fill=category)) +
  geom_bar(stat = "identity", color='black', alpha = 1, width = 0.7) +
  geom_text(aes(label = overlap_sizes, y = overlap_sizes+0.05), 
            size = 5, position = position_dodge(0.9),
            vjust = 0) +
  labs(x = "Intersect sizes",
       y = NULL) +
  theme_classic() +
  scale_fill_manual(values = colors) +
  theme(text = element_text(size= 14, color = "black"),
        legend.position = "none",
        axis.text =element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(hjust = 0),
  ) 
p3
# Overlap matrix
my_overlap_matrix <- str_split(string = my_overlap_sizes$category, pattern = "", simplify = T) %>% 
  as.data.frame() 

colnames(my_overlap_matrix) <- my_names

my_overlap_matrix_tidy <- my_overlap_matrix %>% 
  cbind(category = my_overlap_sizes$category) %>% 
  pivot_longer(cols = !category, names_to = "Set", values_to = "value") %>% 
  full_join(my_overlap_sizes, by = "category") %>% 
  full_join(my_set_sizes, by = "Set")

p4 <- my_overlap_matrix_tidy %>% 
  mutate(category = reorder(category, -overlap_sizes)) %>%  
  mutate(Set = reorder(Set, sizes)) %>%  
  ggplot(aes(x = Set, y = category))+
  geom_tile(aes(fill = Set, alpha = value), color = "grey30", size = 1) +
  scale_fill_manual(values = color, # feel free to use other colors 
                    limits = my_names) +
  scale_alpha_manual(values = c(0.8, 0),  # color the grid for 1, don't color for 0. 
                     limits = c("1", "0")) +
  labs(x = "Sets",  
       y = "Overlap") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(color = "black", size= 14),
        panel.grid = element_blank(),
        axis.text = element_blank()
  )+
  coord_flip()

# put them together
wrap_plots(p1, p2, p4, p3, 
           nrow = 2, 
           ncol = 2,
           heights = c(1, 3), # the more rows in the lower part, the longer it should be
           widths = c(0.3, 0.5),
           guides = "collect") &
  theme(legend.position = "none")


wrap_plots(p3, plot_spacer(), p4, p1, 
           nrow = 2, 
           ncol = 2,
           widths = c(3, 2), # the more rows in the lower part, the longer it should be
           heights = c(0.5, 0.2),
           guides = "collect") &
  theme(legend.position = "none")

ggsave("/Users/suyc/Project/dab/plot/deg/upset.pdf", height = 6, width = 9, bg = "white") 
# this should be a tall & skinny plot 
# I prefer .svg, but you can also save as phd or png 
# I will open up the .svg file and mannually adjust the size until it's good
# check that nothing is cut off from the plot 
# png is for twitter posting 