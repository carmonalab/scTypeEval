# helper to compute not-sample aware (NSA) internal validation metrics

# -silhouette width
# -cLISI
# -ROGUE

nsa_silhouette <- function(scTypeEval,
                           slot = "single-cell"){
   
   
   mat_ident <- scTypeEval@reductions[[slot]]
   emb <- Matrix::t(mat_ident@embeddings)
   ident <- mat_ident@ident[[1]]
   ident_name <- names(mat_ident@ident)
   
   dist <- compute_fast_euclidean(emb)
   
   sil <- compute_silhouette(dist, ident)
   
   # refactor to scTypeEval results
   df_res <- data.frame(celltype = purge_label(names(sil)),
                        measure = as.vector(sil),
                        consistency_metric = "nsa_silhouette",
                        dissimilarity_method = NA,
                        ident = ident_name)
   
   return(df_res)
}




nsa_clisi <- function(scTypeEval,
                      slot = "single-cell",
                      k = 30){
   
   
   mat_ident <- scTypeEval@reductions[[slot]]
   emb <- Matrix::t(mat_ident@embeddings)
   ident <- mat_ident@ident[[1]]
   ident_name <- names(mat_ident@ident)
   celltypes <- levels(ident)
   
   meta_df <- data.frame(celltype = ident)
   lisi_full <- scIntegrationMetrics::compute_lisi(X = emb,
                                                   meta_df,
                                                   label_colnames = "celltype",
                                                   perplexity = k)
   # Normalize LISI: (raw_lisi - 1) / (n_categories - 1)
   n_celltypes <- length(unique(ident[!is.na(ident)]))
   lisi_normalized <- (lisi_full$celltype - 1) / (n_celltypes - 1)
   # Average normalized LISI per cell type
   lisi_scores <- tapply(1 - lisi_normalized, ident, mean, na.rm = TRUE)[celltypes]
   
   # refactor to scTypeEval results
   df_res <- data.frame(celltype = purge_label(names(lisi_scores)),
                        measure = as.vector(lisi_scores),
                        consistency_metric = "nsa_cLISI",
                        dissimilarity_method = NA,
                        ident = ident_name)
   
   return(df_res)
}

nsa_rogue <- function(scTypeEval,
                      slot = "single-cell"){
   
   
   emb <- scTypeEval@reductions[["single-cell"]]@embeddings
   labels <- scTypeEval@reductions[[slot]]@ident[[1]]
   hvg <- scTypeEval@reductions[["single-cell"]]@gene_list
   rogue_counts <- scTypeEval@counts[hvg, colnames(emb)]
   
   rogue_scores <- ROGUE::rogue(rogue_counts,
                                labels = labels,
                                samples = rep("s", length(labels)),
                                platform = "UMI") %>% 
      unlist()
   
   # refactor to scTypeEval results
   df_res <- data.frame(celltype = purge_label(names(rogue_scores)),
                        measure = as.vector(rogue_scores),
                        consistency_metric = "nsa_ROGUE",
                        dissimilarity_method = NA,
                        ident = ident_name)
   
   return(df_res)
}

