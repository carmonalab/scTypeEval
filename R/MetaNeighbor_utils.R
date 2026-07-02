

# Function to keep filtering of min samples and cells but on raw counts
get_filtered_raw_matrix <- function(){
   
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
