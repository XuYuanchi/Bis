# gene pseudotime
re =  readRDS('results/ti/re.rds')
color = c("#F05C3BE5", "#FED439E5", "#46732EE5", "#D5E4A2E5",  "#709AE1E5", "#FD8CC1E5")
p1 = list()
for(i in c(1:11)){
ESC_genes <- c("POU5F1", "NANOG")

ESC_lineage_cds <- re[[i]][[4]][rowData(re[[i]][[4]])$gene_short_name %in% ESC_genes,]

p1[[i]] <- plot_genes_in_pseudotime(ESC_lineage_cds,
                                   color_cells_by="Time",
                                   min_expr=1,
                                   cell_size = 0.5,
                                   label_by_short_name = F,
                                   vertical_jitter = 0,
                                   horizontal_jitter = 0)+
  theme(text = element_text(size=8),legend.position = "none")+
  scale_color_manual(values = color)+
  scale_y_continuous(name="Expression\n(log10)", trans = "log10",
                     breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(.x)))
}
gene_list = c('POU5F1', 'NANOG', 'SOX2', 'DNMT3B', 'NODAL', 'EOMES', 'ID1', 'CDX1', 'MSX2', 'CER1', 'GATA4', 'DKK4',
              'MYCT1', 'POU2AF1','PRDM1', 'CXCR4', 'SOX17', 'YAP1')
Methods = c('Dropout', 'Bis', 'ALRA', 'DCA', 'DeepImpute', 'MAGIC', 'SAVER', 
            'scImpute', 'SCRABBLE', 'scScope', 'scVI')
for(i in c(1:11)){
  #for(j in c(1:length(gene_list))){
  for(j in c(18)){
  ESC_genes <- gene_list[j]
  
  ESC_lineage_cds <- re[[i]][[4]][rowData(re[[i]][[4]])$gene_short_name %in% ESC_genes,]
  
  p2 <- plot_genes_in_pseudotime(ESC_lineage_cds,
                                      color_cells_by="Time",
                                      min_expr=1,
                                      cell_size = 0.5,
                                      label_by_short_name = F,
                                      vertical_jitter = 0,
                                      horizontal_jitter = 0)+
    theme(text = element_text(size=8),legend.position = "none")+
    scale_color_manual(values = color)+
    ggtitle(Methods[i])+
    scale_y_continuous(name="Expression\n(log10)", trans = "log10",
                       breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(.x)))
   ggsave(plot = p2, filename = paste0('plot/ti/', ESC_genes, '_', Methods[i], '.pdf'), height = 2, width = 3)
  }}

