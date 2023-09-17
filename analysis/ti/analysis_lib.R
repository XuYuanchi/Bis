# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="00h"){
  cell_ids <- which(colData(cds)[, "Time"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
# generate trajectories
generate_ICT <- function(cds, cell_metadata, gene_metadata, method){
  cds.raw <- new_cell_data_set(d,
                               cell_metadata = cell_meta,
                               gene_metadata = gene.anno)
  #pre-process the data, using PCA with 100 components
  cds.raw <- preprocess_cds(cds.raw, num_dim = 50)
  
  # reduction the dimension using one of "tSNE", "PCA", "LSI", "Aligned"
  cds.raw <- reduce_dimension(cds.raw, preprocess_method="PCA") # PCA is default, "tSNE", "PCA", "LSI", "Aligned"
  
  cds.raw <- cluster_cells(cds.raw, reduction_method = "UMAP", k=180)
  cds.raw <- learn_graph(cds.raw, use_partition = F, close_loop = F,
                         learn_graph_control = NULL, verbose = FALSE)
  
  p1 <- plot_cells(cds.raw, color_cells_by="Time", group_label_size = 2, cell_size = 0.5,
                       label_cell_groups=F,
                       label_leaves=F,
                       label_branch_points=F,
                       graph_label_size=0)+ggtitle(method)+theme(text = element_text(size=8),legend.position = "none")
  
  cds.raw <- order_cells(cds.raw, root_pr_nodes=get_earliest_principal_node(cds.raw))
  ESC_genes <- c("POU5F1", "NANOG")
  raw.ESC_lineage_cds <- cds.raw[rowData(cds.raw)$gene_short_name %in% ESC_genes,]
  
  p2 <- plot_genes_in_pseudotime(raw.ESC_lineage_cds,
                                     color_cells_by="Time",
                                     min_expr=1,
                                     cell_size = 0.5,
                                     label_by_short_name = F,
                                     vertical_jitter = 0,
                                     horizontal_jitter = 0
  )+theme(text = element_text(size=8),legend.position = "none")+scale_y_continuous(name="Expression\n(log10)", trans = "log10",
                                                                                   breaks = trans_breaks("log10", function(x) 10^x),
                                                                                   labels = trans_format("log10", math_format(.x)))
  
  ESC_genes <- c("CER1", "HNF1B")
  
  raw.ESC_lineage_cds <- cds.raw[rowData(cds.raw)$gene_short_name %in% ESC_genes,]
  
  
  p3 <- plot_genes_in_pseudotime(raw.ESC_lineage_cds,
                                     color_cells_by="Time",
                                     min_expr=0,
                                     cell_size = 0.5,
                                     label_by_short_name = F,
                                     vertical_jitter = 0,
                                     horizontal_jitter = 0)+theme(text = element_text(size=8),legend.position = "none")+scale_y_continuous(name="Expression\n(log10)", trans = "log10",
                                                                                                                                           breaks = trans_breaks("log10", function(x) 10^x),
                                                                                                                                           labels = trans_format("log10", math_format(.x)))
  results <- list()
  results$p1 <- p1
  results$p2 <- p2
  results$p3 <- p3
  return(results)
}