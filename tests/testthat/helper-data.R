# Helper functions for generating test data

#' Generate a synthetic single-cell count matrix with multiple samples and cell types
#'
#' @param n_genes Number of genes (default: 500)
#' @param n_cells_per_sample Number of cells per sample (default: 100)
#' @param n_samples Number of samples (default: 6)
#' @param cell_types Character vector of cell types (default: c("CellA", "CellB", "CellC", "CellD"))
#' @param seed Random seed for reproducibility (default: 42)
#'
#' @return List containing count matrix and metadata
generate_test_data <- function(n_genes = 500,
                               n_cells_per_sample = 100,
                               n_samples = 6,
                               cell_types = c("CellA", "CellB", "CellC", "CellD"),
                               seed = 42) {
  set.seed(seed)
  
  # Generate gene names
  gene_names <- paste0("Gene.", seq_len(n_genes))
  
  # Total cells
  total_cells <- n_cells_per_sample * n_samples
  
  # Generate metadata
  samples <- rep(paste0("Sample", seq_len(n_samples)), each = n_cells_per_sample)
  
  # Assign cell types (evenly distributed across samples)
  n_cell_types <- length(cell_types)
  cells_per_type <- n_cells_per_sample %/% n_cell_types
  remainder <- n_cells_per_sample %% n_cell_types
  
  celltypes <- rep(cell_types, times = c(rep(cells_per_type, n_cell_types - 1), 
                                          cells_per_type + remainder))
  celltypes <- rep(celltypes, times = n_samples)
  
  # Generate synthetic count data
  # Different cell types will have different expression patterns
  count_matrix <- matrix(0, nrow = n_genes, ncol = total_cells)
  rownames(count_matrix) <- gene_names
  colnames(count_matrix) <- paste0("Cell_", seq_len(total_cells))
  
  # Add expression patterns specific to each cell type
  for (i in seq_along(cell_types)) {
    ct <- cell_types[i]
    ct_indices <- which(celltypes == ct)
    
    # Marker genes for this cell type (first 50 genes shifted by cell type)
    marker_start <- ((i - 1) * 50 + 1)
    marker_end <- min(marker_start + 49, n_genes)
    marker_genes <- marker_start:marker_end
    
    # High expression in marker genes
    count_matrix[marker_genes, ct_indices] <- matrix(
      rpois(length(marker_genes) * length(ct_indices), lambda = 50),
      nrow = length(marker_genes)
    )
    
    # Background expression in other genes
    other_genes <- setdiff(seq_len(n_genes), marker_genes)
    count_matrix[other_genes, ct_indices] <- matrix(
      rpois(length(other_genes) * length(ct_indices), lambda = 5),
      nrow = length(other_genes)
    )
  }
  
  # Add sample effects (batch effects)
  for (s in seq_len(n_samples)) {
    sample_cells <- which(samples == paste0("Sample", s))
    # Add a small multiplicative batch effect
    batch_effect <- rnorm(1, mean = 1, sd = 0.1)
    count_matrix[, sample_cells] <- round(count_matrix[, sample_cells] * batch_effect)
  }
  
  # Convert to sparse matrix
  count_matrix <- Matrix::Matrix(count_matrix, sparse = TRUE)
  
  # Create metadata
  # Ensure batch assignment covers all cells properly
  batch_pattern <- rep(c("Batch1", "Batch2"), length.out = n_samples)
  batch <- rep(batch_pattern, each = n_cells_per_sample)
  
  # Ensure condition assignment covers all cells properly
  condition <- rep(c("Control", "Treatment"), length.out = total_cells)
  
  metadata <- data.frame(
    celltype = celltypes,
    sample = samples,
    batch = batch,
    condition = condition,
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- colnames(count_matrix)
  
  return(list(counts = count_matrix, metadata = metadata))
}


#' Generate a small test dataset for quick tests
#'
#' @return List containing count matrix and metadata
generate_small_test_data <- function() {
  generate_test_data(
    n_genes = 200,
    n_cells_per_sample = 40,
    n_samples = 5,
    cell_types = c("CellA", "CellB", "CellC"),
    seed = 123
  )
}


#' Create a minimal Seurat object from test data (if Seurat is available)
#'
#' @param test_data List with counts and metadata
#' @return Seurat object or NULL if Seurat is not available
create_test_seurat <- function(test_data = NULL) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    return(NULL)
  }
  
  if (is.null(test_data)) {
    test_data <- generate_small_test_data()
  }
  
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = test_data$counts,
    meta.data = test_data$metadata
  )
  
  return(seurat_obj)
}


#' Create a minimal SingleCellExperiment object from test data
#'
#' @param test_data List with counts and metadata
#' @return SingleCellExperiment object or NULL if not available
create_test_sce <- function(test_data = NULL) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    return(NULL)
  }
  
  if (is.null(test_data)) {
    test_data <- generate_small_test_data()
  }
  
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = test_data$counts),
    colData = test_data$metadata
  )
  
  return(sce)
}


#' Create a complete scTypeEval object for testing
#'
#' @param small Logical, whether to use small dataset (default: TRUE)
#' @return scTypeEval object
create_test_scTypeEval <- function(small = TRUE) {
  if (small) {
    test_data <- generate_small_test_data()
  } else {
    test_data <- generate_test_data()
  }
  
  sceval <- create.scTypeEval(
    matrix = test_data$counts,
    metadata = test_data$metadata
  )
  
  return(sceval)
}


#' Create a fully processed scTypeEval object for downstream testing
#'
#' @param small Logical, whether to use small dataset (default: TRUE)
#' @return Processed scTypeEval object with HVGs computed
create_processed_scTypeEval <- function(small = TRUE) {
  sceval <- create_test_scTypeEval(small = small)
  
  # Run processing
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    min.samples = 3,
    min.cells = 5,
    verbose = FALSE
  )
  
  # Add HVG genes (required for downstream analyses)
  # Use "basic" method to avoid parallelization issues in tests
  sceval <- Run.HVG(
    sceval,
    var.method = "scran",
    verbose = FALSE
  )
  
  return(sceval)
}
