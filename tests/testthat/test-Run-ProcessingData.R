test_that("Run.ProcessingData processes single-cell data", {
  sceval <- create_test_scTypeEval()
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    aggregation = "single-cell",
    verbose = FALSE
  )
  
  expect_true("single-cell" %in% names(sceval@data))
  expect_s4_class(sceval@data[["single-cell"]], "DataAssay")
  expect_s4_class(sceval@data[["single-cell"]]@matrix, "dgCMatrix")
})


test_that("Run.ProcessingData processes pseudobulk data", {
  sceval <- create_test_scTypeEval()
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    aggregation = "pseudobulk",
    verbose = FALSE
  )
  
  expect_true("pseudobulk" %in% names(sceval@data))
  expect_s4_class(sceval@data[["pseudobulk"]], "DataAssay")
  expect_s4_class(sceval@data[["pseudobulk"]]@matrix, "dgCMatrix")
})


test_that("Run.ProcessingData processes both aggregation types", {
  sceval <- create_test_scTypeEval()
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    aggregation = c("single-cell", "pseudobulk"),
    verbose = FALSE
  )
  
  expect_true(all(c("single-cell", "pseudobulk") %in% names(sceval@data)))
})


test_that("Run.ProcessingData uses active.ident when ident is NULL", {
  sceval <- create_test_scTypeEval()
  sceval <- set.activeIdent(sceval, ident = "celltype")
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = NULL,
    sample = "sample",
    verbose = FALSE
  )
  
  expect_true("single-cell" %in% names(sceval@data))
})


test_that("Run.ProcessingData respects min.samples parameter", {
  test_data <- generate_test_data(n_samples = 7)
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    min.samples = 5,
    min.cells = 5,
    verbose = FALSE
  )
  
  # Should filter out cell types with fewer than 2 samples
  cell_types <- unique(sceval@data[["single-cell"]]@ident[[1]])
  expect_true(length(cell_types) > 0)
})


test_that("Run.ProcessingData warns when dataset has fewer than 5 samples", {
  # Create dataset with only 3 samples - should trigger warning
  test_data <- generate_test_data(n_samples = 3)
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  
  # Should produce warning about total number of samples < 5
  expect_warning(
    sceval <- Run.ProcessingData(
      sceval,
      ident = "celltype",
      sample = "sample",
      min.samples = 2,
      min.cells = 5,
      verbose = FALSE
    ),
    "Only 3 samples detected.*For inter-sample comparison 5 or more samples is recommended"
  )
})


test_that("Run.ProcessingData works without warning when dataset has 5 or more samples", {
  # Create dataset with 6 samples - should not trigger warning
  test_data <- generate_test_data(n_samples = 6)
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  
  # Should process successfully without the sample count warning
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    min.samples = 3,
    min.cells = 5,
    verbose = FALSE
  )
  
  expect_true("single-cell" %in% names(sceval@data))
})


test_that("Run.ProcessingData respects min.cells parameter", {
  test_data <- generate_small_test_data()
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    min.samples = 3,
    min.cells = 5,
    verbose = FALSE
  )
  
  # Should process successfully with reasonable thresholds
  expect_true("single-cell" %in% names(sceval@data))
})


test_that("Run.ProcessingData accepts normalization.method parameter", {
  sceval <- create_test_scTypeEval()
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    normalization.method = "Log1p",
    verbose = FALSE
  )
  
  expect_s4_class(sceval@data[["single-cell"]]@matrix, "dgCMatrix")
})


# Tests for different normalization methods ---------------------------------

test_that("Run.ProcessingData works with Log1p normalization", {
  sceval <- create_test_scTypeEval()
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    normalization.method = "Log1p",
    aggregation = "single-cell",
    verbose = FALSE
  )
  
  expect_true("single-cell" %in% names(sceval@data))
  expect_s4_class(sceval@data[["single-cell"]]@matrix, "dgCMatrix")
  
  # Normalized values should be non-negative
  expect_true(all(sceval@data[["single-cell"]]@matrix@x >= 0))
})


test_that("Run.ProcessingData works with CLR normalization", {
  sceval <- create_test_scTypeEval()
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    normalization.method = "CLR",
    aggregation = "single-cell",
    verbose = FALSE
  )
  
  expect_true("single-cell" %in% names(sceval@data))
  expect_s4_class(sceval@data[["single-cell"]]@matrix, "dgCMatrix")
})


test_that("Run.ProcessingData works with pearson normalization", {
  sceval <- create_test_scTypeEval()
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    normalization.method = "pearson",
    aggregation = "single-cell",
    verbose = FALSE
  )
  
  expect_true("single-cell" %in% names(sceval@data))
  expect_s4_class(sceval@data[["single-cell"]]@matrix, "dgCMatrix")
})


test_that("Run.ProcessingData normalization methods produce different results", {
  sceval <- create_test_scTypeEval()
  
  # Log1p normalization
  sceval_log <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    normalization.method = "Log1p",
    verbose = FALSE
  )
  
  # CLR normalization
  sceval_clr <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    normalization.method = "CLR",
    verbose = FALSE
  )
  
  # Pearson normalization
  sceval_pears <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    normalization.method = "pearson",
    verbose = FALSE
  )
  
  # Get matrices
  mat_log <- sceval_log@data[["single-cell"]]@matrix
  mat_clr <- sceval_clr@data[["single-cell"]]@matrix
  mat_pears <- sceval_pears@data[["single-cell"]]@matrix
  
  # All should be dgCMatrix
  expect_s4_class(mat_log, "dgCMatrix")
  expect_s4_class(mat_clr, "dgCMatrix")
  expect_s4_class(mat_pears, "dgCMatrix")
  
  # Results should differ between methods
  expect_false(identical(mat_log, mat_clr))
  expect_false(identical(mat_log, mat_pears))
  expect_false(identical(mat_clr, mat_pears))
})


test_that("Run.ProcessingData Log1p normalization with pseudobulk", {
  sceval <- create_test_scTypeEval()
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    normalization.method = "Log1p",
    aggregation = "pseudobulk",
    verbose = FALSE
  )
  
  expect_true("pseudobulk" %in% names(sceval@data))
  expect_s4_class(sceval@data[["pseudobulk"]]@matrix, "dgCMatrix")
  expect_true(all(sceval@data[["pseudobulk"]]@matrix@x >= 0))
})


test_that("Run.ProcessingData CLR normalization with pseudobulk", {
  sceval <- create_test_scTypeEval()
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    normalization.method = "CLR",
    aggregation = "pseudobulk",
    verbose = FALSE
  )
  
  expect_true("pseudobulk" %in% names(sceval@data))
  expect_s4_class(sceval@data[["pseudobulk"]]@matrix, "dgCMatrix")
})


test_that("Run.ProcessingData pearson normalization with pseudobulk", {
  sceval <- create_test_scTypeEval()
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    normalization.method = "pearson",
    aggregation = "pseudobulk",
    verbose = FALSE
  )
  
  expect_true("pseudobulk" %in% names(sceval@data))
  expect_s4_class(sceval@data[["pseudobulk"]]@matrix, "dgCMatrix")
})


test_that("Run.ProcessingData normalization methods with both aggregations", {
  sceval <- create_test_scTypeEval()
  
  # Test each normalization method with both aggregations
  for (norm_method in c("Log1p", "CLR", "pearson")) {
    sceval_test <- Run.ProcessingData(
      sceval,
      ident = "celltype",
      sample = "sample",
      normalization.method = norm_method,
      aggregation = c("single-cell", "pseudobulk"),
      verbose = FALSE
    )
    
    expect_true(all(c("single-cell", "pseudobulk") %in% names(sceval_test@data)),
                info = paste("Both aggregations should exist for", norm_method))
    expect_s4_class(sceval_test@data[["single-cell"]]@matrix, "dgCMatrix")
    expect_s4_class(sceval_test@data[["pseudobulk"]]@matrix, "dgCMatrix")
  }
})


test_that("Run.ProcessingData invalid normalization method errors gracefully", {
  sceval <- create_test_scTypeEval()
  
  expect_error(
    Run.ProcessingData(
      sceval,
      ident = "celltype",
      sample = "sample",
      normalization.method = "invalid_method",
      verbose = FALSE
    ),
    "not a supported normalization method"
  )
})


test_that("Run.ProcessingData normalization methods preserve matrix dimensions", {
  sceval <- create_test_scTypeEval()
  
  # Get dimensions before normalization
  n_genes <- nrow(sceval@counts)
  n_cells <- ncol(sceval@counts)
  
  for (norm_method in c("Log1p", "CLR", "pearson")) {
    sceval_test <- Run.ProcessingData(
      sceval,
      ident = "celltype",
      sample = "sample",
      normalization.method = norm_method,
      aggregation = "single-cell",
      min.samples = 2,
      min.cells = 5,
      verbose = FALSE
    )
    
    # Number of genes should be preserved
    expect_equal(nrow(sceval_test@data[["single-cell"]]@matrix), n_genes)
    
    # Number of cells might be reduced due to filtering, but should be positive
    expect_true(ncol(sceval_test@data[["single-cell"]]@matrix) > 0,
               info = paste("Should have cells after filtering for", norm_method))
  }
})


test_that("Run.ProcessingData stores correct metadata in DataAssay", {
  sceval <- create_test_scTypeEval()
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    verbose = FALSE
  )
  
  data_assay <- sceval@data[["single-cell"]]
  expect_equal(data_assay@aggregation, "single-cell")
  expect_true(is.factor(data_assay@group))
  expect_true(is.factor(data_assay@sample))
  expect_true(is.list(data_assay@ident))
  expect_true("celltype" %in% names(data_assay@ident))
})


test_that("Run.ProcessingData handles verbose parameter", {
  sceval <- create_test_scTypeEval()
  
  # verbose = TRUE should produce messages
  expect_message(
    sceval <- Run.ProcessingData(
      sceval,
      ident = "celltype",
      sample = "sample",
      verbose = TRUE
    ),
    "Processing data"
  )
  
  # verbose = FALSE should suppress messages
  sceval2 <- create_test_scTypeEval()
  expect_silent(
    sceval2 <- Run.ProcessingData(
      sceval2,
      ident = "celltype",
      sample = "sample",
      verbose = FALSE
    )
  )
})


test_that("Run.ProcessingData errors on invalid ident", {
  sceval <- create_test_scTypeEval()
  
  expect_error(
    Run.ProcessingData(
      sceval,
      ident = "nonexistent_column",
      sample = "sample",
      verbose = FALSE
    )
  )
})


test_that("Run.ProcessingData errors on invalid sample", {
  sceval <- create_test_scTypeEval()
  
  expect_error(
    Run.ProcessingData(
      sceval,
      ident = "celltype",
      sample = "nonexistent_column",
      verbose = FALSE
    )
  )
})


test_that("Run.ProcessingData errors on invalid aggregation type", {
  sceval <- create_test_scTypeEval()
  
  expect_error(
    Run.ProcessingData(
      sceval,
      ident = "celltype",
      sample = "sample",
      aggregation = "invalid_type",
      verbose = FALSE
    ),
    "Invalid aggregation type"
  )
})


test_that("Run.ProcessingData handles multiple samples correctly", {
  test_data <- generate_test_data(n_samples = 6)
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    min.samples = 5,
    min.cells = 10,
    verbose = FALSE
  )
  
  # Should retain cell types present in at least 5 samples
  expect_true("single-cell" %in% names(sceval@data))
  expect_true(ncol(sceval@data[["single-cell"]]@matrix) > 0)
})


test_that("Run.ProcessingData creates proper pseudobulk aggregation", {
  test_data <- generate_test_data(n_samples = 6)
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    aggregation = "pseudobulk",
    min.samples = 3,
    min.cells = 10,
    verbose = FALSE
  )
  
  pb_assay <- sceval@data[["pseudobulk"]]
  
  # Pseudobulk should have one column per sample-celltype combination
  expect_true(ncol(pb_assay@matrix) > 0)
  expect_equal(length(pb_assay@sample), ncol(pb_assay@matrix))
  expect_equal(length(pb_assay@ident[[1]]), ncol(pb_assay@matrix))
})


test_that("Run.ProcessingData errors when min.samples is too high", {
  test_data <- generate_test_data(n_samples = 5)
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  
  # Require more samples than available - should result in error or empty data
  expect_error(
    Run.ProcessingData(
      sceval,
      ident = "celltype",
      sample = "sample",
      min.samples = 10,  # More than the 3 samples available
      min.cells = 5,
      verbose = FALSE
    )
  )
})


test_that("Run.ProcessingData errors when min.cells is too high", {
  test_data <- generate_test_data(n_samples = 5)
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  
  # Require more cells per sample-celltype than available
  expect_error(
    Run.ProcessingData(
      sceval,
      ident = "celltype",
      sample = "sample",
      min.samples = 5,
      min.cells = 1000,  # More cells than available per group
      verbose = FALSE
    )
  )
})


test_that("Run.ProcessingData filtering reduces cell types correctly", {
  # Create data with 6 samples and 4 cell types
  test_data <- generate_test_data(
    n_samples = 6,
    n_cells_per_sample = 100,
    cell_types = c("CellA", "CellB", "CellC", "CellD")
  )
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  
  # Count original cell types
  original_celltypes <- unique(test_data$metadata$celltype)
  
  # Process with lenient thresholds - should keep all cell types
  sceval_lenient <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    min.samples = 3,
    min.cells = 10,
    verbose = FALSE
  )
  
  retained_celltypes <- unique(sceval_lenient@data[["single-cell"]]@ident[[1]])
  
  # All cell types should be retained with lenient thresholds
  expect_equal(length(retained_celltypes), length(original_celltypes))
  expect_true(all(sort(as.character(retained_celltypes)) == sort(original_celltypes)))
})


test_that("Run.ProcessingData filtering verifies correct number of samples per cell type", {
  # Create data with known structure
  test_data <- generate_test_data(n_samples = 6, n_cells_per_sample = 80)
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    min.samples = 4,
    min.cells = 10,
    verbose = FALSE
  )
  
  # Check that each retained cell type appears in at least min.samples
  sc_data <- sceval@data[["single-cell"]]
  
  # Extract sample-celltype combinations
  for (ct in unique(sc_data@ident[[1]])) {
    # Get all groups for this cell type
    ct_groups <- sc_data@group[grepl(paste0("_", ct, "$"), sc_data@group)]
    # Extract unique samples for this cell type
    ct_samples <- unique(sapply(as.character(ct_groups), function(x) strsplit(x, "_")[[1]][1]))
    
    # Each cell type should be in at least min.samples (4)
    expect_true(length(ct_samples) >= 4,
                info = paste("Cell type", ct, "has", length(ct_samples), "samples"))
  }
})


test_that("Run.ProcessingData filtering verifies correct number of cells per group", {
  # Create data with known structure
  test_data <- generate_test_data(n_samples = 6, n_cells_per_sample = 100)
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    min.samples = 3,
    min.cells = 15,
    verbose = FALSE
  )
  
  # Check that each sample-celltype combination has at least min.cells
  sc_data <- sceval@data[["single-cell"]]
  
  # Count cells per group
  group_counts <- table(sc_data@group)
  
  # Each group should have at least min.cells (15)
  expect_true(all(group_counts >= 15),
              info = paste("Some groups have fewer than 15 cells:",
                          paste(names(group_counts[group_counts < 15]), collapse = ", ")))
})


test_that("Run.ProcessingData pseudobulk filtering verifies sample-celltype combinations", {
  test_data <- generate_test_data(n_samples = 6, n_cells_per_sample = 80)
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    aggregation = "pseudobulk",
    min.samples = 4,
    min.cells = 10,
    verbose = FALSE
  )
  
  pb_data <- sceval@data[["pseudobulk"]]
  
  # Each cell type should appear in at least min.samples (4)
  for (ct in unique(pb_data@ident[[1]])) {
    ct_samples <- pb_data@sample[pb_data@ident[[1]] == ct]
    n_samples <- length(unique(ct_samples))
    
    expect_true(n_samples >= 4,
                info = paste("Cell type", ct, "appears in", n_samples, "samples"))
  }
})


test_that("Run.ProcessingData correctly removes cell types below threshold", {
  # Create custom data where we know one cell type will be filtered
  test_data <- generate_test_data(
    n_samples = 6,
    n_cells_per_sample = 100,
    cell_types = c("CellA", "CellB", "CellC", "CellD")
  )
  
  # Manually remove CellD from some samples to make it appear in fewer samples
  metadata <- test_data$metadata
  # Keep CellD only in first 2 samples
  to_remove <- metadata$celltype == "CellD" & 
               !metadata$sample %in% c("Sample1", "Sample2")
  
  test_data$counts <- test_data$counts[, !to_remove]
  test_data$metadata <- test_data$metadata[!to_remove, ]
  
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    min.samples = 4,  # CellD only in 2 samples, should be filtered
    min.cells = 5,
    verbose = FALSE
  )
  
  retained_celltypes <- unique(sceval@data[["single-cell"]]@ident[[1]])
  
  # CellD should be removed
  expect_false("CellD" %in% retained_celltypes,
               info = "CellD should be filtered out (only in 2 samples, threshold is 4)")
  
  # Other cell types should remain
  expect_true(all(c("CellA", "CellB", "CellC") %in% retained_celltypes))
})
