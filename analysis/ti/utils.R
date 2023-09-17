# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="zygote"){
  cell_ids <- which(colData(cds)[, "Time"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}


process <- function(d, cell_meta, gene.anno, methods_name, cell_label, time_bin){
  cds.raw <- new_cell_data_set(d,
                               cell_metadata = cell_meta,
                               gene_metadata = gene.anno)
  #pre-process the data, using PCA with 100 components
  cds.raw <- preprocess_cds(cds.raw, num_dim = 50)
  
  # reduction the dimension using one of "tSNE", "PCA", "LSI", "Aligned"
  cds.raw <- reduce_dimension(cds.raw, preprocess_method="PCA") # PCA is default, "tSNE", "PCA", "LSI", "Aligned"
  
  cds.raw <- cluster_cells(cds.raw, reduction_method = "UMAP")
  cds.raw <- learn_graph(cds.raw, use_partition = F, close_loop = F,
                         learn_graph_control = NULL, verbose = FALSE)
  
  raw.p1 <- plot_cells(cds.raw, color_cells_by="Time", group_label_size = 6, cell_size = 2,
                       label_cell_groups=F,
                       label_leaves=F,
                       label_branch_points=F,
                       graph_label_size=4)+ggtitle(methods_name)+theme(text = element_text(size=8),legend.position = "bottom")
  
  cds.raw <- order_cells(cds.raw, root_pr_nodes=get_earliest_principal_node(cds.raw, time_bin))
  
  monocle3_pseudotime = cds.raw@principal_graph_aux[['UMAP']]$pseudotime
  monocle3_pseudotime[which(!is.finite(monocle3_pseudotime))] = 0
  cor.kendall = cor(monocle3_pseudotime, as.numeric(cell_label), method = "kendall", use = "complete.obs")
  
  # cor.kendall = cor(cds.raw@phenoData@data$Pseudotime, as.numeric(cell_label), method = "kendall", use = "complete.obs")
  
  # cor.kendall = cor(cds.raw@principal_graph_aux@listData$UMAP$pseudotime, as.numeric(cds.raw@colData@listData$Time), 
  #                   method = "kendall", use = "na.or.complete")
  
  lpsorder2 = data.frame(sample_name = cell_meta$Time, State= cds.raw@colData@listData$Time,
                         Pseudotime = cds.raw@principal_graph_aux@listData$UMAP$pseudotime, rank = rank(cds.raw@principal_graph_aux@listData$UMAP$pseudotime))
  
  lpsorder_rank = dplyr::arrange(lpsorder2, rank)
  
  lpsorder_rank$Pseudotime = lpsorder_rank$rank
  
  lpsorder_rank = lpsorder_rank[-4]
  
  lpsorder_rank[1] <- lapply(lpsorder_rank[1], as.character)
  
  subpopulation <- data.frame(cell = cell_meta$Time, sub = as.numeric(cell_label)-1)
  
  POS <- TSCAN::orderscore(subpopulation, lpsorder_rank['sample_name'])[[1]]
  
  Results = list()
  Results[[1]] = raw.p1
  Results[[2]] = cor.kendall
  Results[[3]] = POS
  Results[[4]] = cds.raw
  return(Results)
}

add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}