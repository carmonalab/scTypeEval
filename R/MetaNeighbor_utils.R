

# Function to keep filtering of min samples and cells but on raw counts
get_filtered_raw_matrix <- function(scTypeEval){
   
   if(is.null(sceval@data$`single-cell`)){
      stop("No normalization slot found. Please run before `run_processing_data()`.\n")
   }
   
   k <- colnames(sceval@data$`single-cell`@matrix)
   sub <- list(counts = sceval@counts[,k],
               metadata = sceval@metadata[k,])
   return(sub)
}
   
   
get_sce <- function(counts,
                    cell_metadata = NULL, # cell metadata
                    gene_metadata = NULL, # gene metadata
                    coldata_cols = NULL # which cols to keep
                    ) {
   
   if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Package 'SummarizedExperiment' is required.")
   }
   
   # Basic checks
   if (is.null(rownames(counts)) || is.null(colnames(counts))) {
      stop("Counts matrix must have rownames (genes) and colnames (cells).")
   }
   
   # Build colData (cell metadata)
   if (!is.null(cell_metadata)) {
      cell_metadata <- as.data.frame(cell_metadata)
      
      # ensure matching order with counts columns
      cell_metadata <- cell_metadata[colnames(counts), , drop = FALSE]
      
      if (!is.null(coldata_cols)) {
         cell_metadata <- cell_metadata[, coldata_cols, drop = FALSE]
      }
   } else {
      cell_metadata <- data.frame(row.names = colnames(counts))
   }
   
   # Build rowData (gene metadata)
   if (!is.null(gene_metadata)) {
      gene_metadata <- as.data.frame(gene_metadata)
      
      # ensure matching order with counts rows
      gene_metadata <- gene_metadata[rownames(counts), , drop = FALSE]
   } else {
      gene_metadata <- data.frame(row.names = rownames(counts))
   }
   
   se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts),
      colData = cell_metadata,
      rowData = gene_metadata
   )
   
   return(se)
}



# Utility to get variable genes according to MetaNeighbor workflow
MN_variableGenes <- function(mat,
                             exp_labels){
   sce <- get_sce(counts = mat)
   
   var_genes = MetaNeighbor::variableGenes(dat = sce,
                                           exp_labels = exp_labels)
   
   return(var_genes)
}

# returns directly value of consistency per cell type
MN_supervised <- function(scTypeEval,
                          gene_list = NULL,
                          verbose = TRUE)
{
   
   # filter cells that only passed the processing
   filt_list <- get_filtered_raw_matrix(scTypeEval)
   # convert to sce
   sce <- get_sce(counts = filt_list$counts,
                  cell_metadata = filt_list$metadata)
   # get sample and ident from processed metadata
   exp_labs <- scTypeEval@data$`single-cell`@sample
   ident <- names(scTypeEval@data$`single-cell`@ident)
   # Metaneighbor required conversion of celltypes into a oneshot dataframe
   labs <- filt_list$metadata[[ident]]
   # Create one-hot encoding
   one_hot <- model.matrix(~ labs - 1)
   # Clean column names
   colnames(one_hot) <- sub("^labs", "", colnames(one_hot))
   # Preserve row names
   rownames(one_hot) <- rownames(filt_list$metadata)
   
   geneset <- check_genelist(scTypeEval = scTypeEval,
                             gene_list = gene_list,
                             verbose = verbose)
   
   res <- MetaNeighbor::MetaNeighbor(dat = sce,
                                     i = 1,
                                     experiment_labels = exp_labs,
                                     celltype_labels = one_hot,
                                     genesets = list(GeneSet = geneset),
                                     bplot = F,
                                     fast_version = T,
                                     node_degree_normalization = T,
                                     batch_size = 10,
                                     detailed_results = F)
   # refactor to scTypeEval results
   df_res <- data.frame(celltype = purge_label(colnames(res)),
                        measure = as.vector(res),
                        consistency_metric = "MetaNeighbor_Supervised",
                        dissimilarity_method = NA,
                        ident = ident)
   
   return(df_res)
   
}


MN_unsupervised <- function(scTypeEval,
                            group,
                            gene_list = NULL,
                            method = "score",
                            verbose = TRUE){
   # filter cells that only passed the processing
   filt_list <- get_filtered_raw_matrix(scTypeEval)
   # convert to sce
   sce <- get_sce(counts = filt_list$counts,
                  cell_metadata = filt_list$metadata)
   # get sample and ident from processed metadata
   exp_labs <- scTypeEval@data$`single-cell`@sample
   ident <- scTypeEval@data$`single-cell`@ident[[1]]
   
   geneset <- check_genelist(scTypeEval = scTypeEval,
                             gene_list = gene_list,
                             verbose = verbose)
   
   res_score <- MetaNeighbor::MetaNeighborUS(var_genes = geneset,
                                       dat = sce,
                                       i = 1,
                                       study_id = exp_labs,
                                       cell_type = ident,
                                       trained_model = NULL,
                                       fast_version = TRUE,
                                       node_degree_normalization = TRUE,
                                       one_vs_best = FALSE,
                                       symmetric_output = TRUE)
   
   # use MetaNeighbor function to find the reciprocal top hits in each target study
   # but keeping self comparisons
   top_hits <- topHitsByStudy(auroc = res_score,
                              threshold = 0,
                              n_digits = 3,
                              collapse_duplicates = T) |>
      mutate(match_numeric = ifelse(Match_type == "Reciprocal_top_hit",
                                    0, 1))
   # generate new matrix to populate with reciprocal_top_hit (0) or not (1)
   res_match <- matrix(
      1,
      nrow = nrow(res_score),
      ncol = ncol(res_score),
      dimnames = list(rownames(res_score), colnames(res_score))
   )
   
   # Fill in symmetric matrix
   for (a in seq_len(nrow(top_hits))) {
      tmp <- top_hits[a,]
      i <- tmp[[1]]
      j <- tmp[[2]]
      d <- tmp[["match_numeric"]]
      res_match[i, j] <- res_match[j, i] <- as.numeric(d)
   }
   
   rownames(res_match) <- gsub("[|]", "_", rownames(res_match))
   colnames(res_match) <- gsub("[|]", "_", colnames(res_match))
   
   res <- switch (method,
                  "score" = 1 - res_score, # invert scale
                  "match" = res_match,
                  stop(method, "is not a supported method for MetaNeighbor.\n")
   )
   
   # adapt row/colnames for compatability with scTypeEval
   rownames(res) <- gsub("[|]", "_", rownames(res))
   colnames(res) <- gsub("[|]", "_", colnames(res))
   
   # adjust ordering
   res <- res[group, group]
   
   return(as.dist(res))
   
}

# reuse function from MetaNeighbor to not filter self comparisons
# https://github.com/gillislab/MetaNeighbor/blob/master/R/topHits.R
#87e4e7a36bab3b6cb38c58d6bd0117c93dc69e9c

topHitsByStudy = function(auroc, threshold = 0.9, n_digits = 2, collapse_duplicates = TRUE) {
   `%>%` <- dplyr::`%>%`
   #could be sped up by finding same study AUROC's when in matrix form and masking them to NA (as in the old topHits function)
   result <- tibble::as_tibble(auroc, rownames = "ref_cell_type") %>%
      tidyr::pivot_longer(cols = -ref_cell_type,
                          names_to = "target_cell_type",
                          values_to = "auroc") %>%
      dplyr::mutate(ref_study = getStudyId(ref_cell_type),
                    target_study = getStudyId(target_cell_type)) %>%
      # CHANGE: keep self - remove filter
      #dplyr::filter(ref_study != target_study) %>%
      dplyr::group_by(ref_cell_type, target_study) %>%
      dplyr::filter(auroc == max(auroc, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-ref_study, -target_study) %>%
      dplyr::mutate(is_reciprocal = is_reciprocal_top_hit(.)) %>%
      dplyr::filter(auroc >= threshold) 
   
   if (collapse_duplicates) {
      result <- result %>%
         dplyr::group_by(ref_cell_type, target_cell_type) %>%
         dplyr::mutate(pair_id = paste(sort(c(ref_cell_type, target_cell_type)),
                                       collapse = "")) %>%
         dplyr::group_by(pair_id) %>%
         dplyr::summarize(ref_cell_type = dplyr::first(ref_cell_type),
                          target_cell_type = dplyr::first(target_cell_type),
                          auroc = mean(auroc),
                          is_reciprocal = dplyr::first(is_reciprocal)) %>%
         dplyr::ungroup() %>%
         dplyr::select(-pair_id)
   }
   
   # final formatting
   result <- result %>%
      dplyr::arrange(desc(auroc)) %>%
      dplyr::mutate(auroc = round(auroc, n_digits)) %>%    
      dplyr::mutate(Match_type = ifelse(is_reciprocal,
                                        "Reciprocal_top_hit",
                                        paste0("Above_", threshold))) %>%
      dplyr::select("Study_ID|Celltype_1" = ref_cell_type,
                    "Study_ID|Celltype_2" = target_cell_type,
                    "AUROC" = auroc,
                    Match_type)
   return(result)
}



is_reciprocal_top_hit = function(best_hits) {
   `%>%` <- dplyr::`%>%`
   best_hits <- best_hits %>%
      dplyr::select(-auroc)
   reverse_hits <- best_hits %>%
      dplyr::select(target_cell_type = ref_cell_type,
                    reciprocal_cell_type = target_cell_type)
   reciprocal_best_hits <- dplyr::inner_join(best_hits, reverse_hits,
                                             by = "target_cell_type") %>%
      dplyr::filter(ref_cell_type == reciprocal_cell_type) %>%
      dplyr::select(-reciprocal_cell_type) %>%
      tibble::add_column(is_reciprocal = TRUE)
   result <- dplyr::left_join(best_hits, reciprocal_best_hits,
                              by = c("ref_cell_type", "target_cell_type")) %>%
      tidyr::replace_na(replace = list(is_reciprocal = FALSE)) %>%
      dplyr::pull(is_reciprocal)
   return(result)
}
