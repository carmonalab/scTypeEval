

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
   exp_labs <- sceval@data$`single-cell`@sample
   ident <- names(sceval@data$`single-cell`@ident)
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
