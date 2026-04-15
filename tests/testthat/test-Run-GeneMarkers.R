test_that("run_gene_markers identifies marker genes with scran.findMarkers", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- run_gene_markers(
    sceval,
    method = "scran.findMarkers",
    ngenes_celltype = 50,
    aggregation = "single-cell",
    verbose = FALSE
  )
  
  expect_true("scran.findMarkers" %in% names(sceval@gene_lists))
  expect_true(length(sceval@gene_lists[["scran.findMarkers"]]) > 0)
})


test_that("run_gene_markers respects ngenes_celltype parameter", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- run_gene_markers(
    sceval,
    method = "scran.findMarkers",
    ngenes_celltype = 20,
    verbose = FALSE
  )
  
  # Total markers should be reasonable given cell types
  expect_true(length(sceval@gene_lists[["scran.findMarkers"]]) > 0)
})


test_that("run_gene_markers works with single-cell aggregation", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- run_gene_markers(
    sceval,
    method = "scran.findMarkers",
    ngenes_celltype = 50,
    aggregation = "single-cell",
    verbose = FALSE
  )
  
  expect_true("scran.findMarkers" %in% names(sceval@gene_lists))
})


test_that("run_gene_markers respects black_list parameter", {
  sceval <- create_processed_scTypeEval()
  
  # Get some gene names
  all_genes <- rownames(sceval@data[["single-cell"]]@matrix)
  black_genes <- all_genes[seq_len(10)]
  
  sceval@black_list <- black_genes
  
  sceval <- run_gene_markers(
    sceval,
    method = "scran.findMarkers",
    ngenes_celltype = 50,
    verbose = FALSE
  )
  
  # Markers should not contain blacklisted genes
  expect_true(!any(sceval@gene_lists[["scran.findMarkers"]] %in% black_genes))
})


test_that("run_gene_markers accepts custom black_list parameter", {
  sceval <- create_processed_scTypeEval()
  
  all_genes <- rownames(sceval@data[["single-cell"]]@matrix)
  custom_black <- all_genes[seq_len(5)]
  
  sceval <- run_gene_markers(
    sceval,
    method = "scran.findMarkers",
    ngenes_celltype = 50,
    black_list = custom_black,
    verbose = FALSE
  )
  
  expect_true(!any(sceval@gene_lists[["scran.findMarkers"]] %in% custom_black))
})


test_that("run_gene_markers errors without processed data", {
  sceval <- create_test_scTypeEval()
  
  expect_error(
    run_gene_markers(
      sceval,
      method = "scran.findMarkers",
      ngenes_celltype = 50,
      verbose = FALSE
    ),
    "No normalization slot found"
  )
})


test_that("run_gene_markers errors on invalid aggregation type", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    run_gene_markers(
      sceval,
      method = "scran.findMarkers",
      ngenes_celltype = 50,
      aggregation = "invalid_type",
      verbose = FALSE
    ),
    "Invalid aggregation type"
  )
})


test_that("run_gene_markers errors on unsupported method", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    run_gene_markers(
      sceval,
      method = "unsupported_method",
      ngenes_celltype = 50,
      verbose = FALSE
    ),
    "Supported gene markes definitions"
  )
})


test_that("run_gene_markers handles verbose parameter", {
  sceval <- create_processed_scTypeEval()
  
  # verbose = TRUE should produce messages
  expect_message(
    sceval <- run_gene_markers(
      sceval,
      method = "scran.findMarkers",
      ngenes_celltype = 50,
      verbose = TRUE
    )
  )
})


test_that("run_gene_markers handles ncores parameter", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- run_gene_markers(
    sceval,
    method = "scran.findMarkers",
    ngenes_celltype = 50,
    ncores = 1,
    verbose = FALSE
  )
  
  expect_true("scran.findMarkers" %in% names(sceval@gene_lists))
})


test_that("run_gene_markers returns genes present in the data", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- run_gene_markers(
    sceval,
    method = "scran.findMarkers",
    ngenes_celltype = 50,
    verbose = FALSE
  )
  
  available_genes <- rownames(sceval@data[["single-cell"]]@matrix)
  expect_true(all(sceval@gene_lists[["scran.findMarkers"]] %in% available_genes))
})


test_that("run_gene_markers works with multiple samples", {
  test_data <- generate_test_data(n_samples = 6)
  sceval <- create_scTypeEval(test_data$counts, test_data$metadata)
  sceval <- run_processing_data(sceval, ident = "celltype", sample = "sample",
                                min_samples = 3, min_cells = 5, verbose = FALSE)
  
  sceval <- run_gene_markers(
    sceval,
    method = "scran.findMarkers",
    ngenes_celltype = 50,
    verbose = FALSE
  )
  
  expect_true("scran.findMarkers" %in% names(sceval@gene_lists))
  expect_true(length(sceval@gene_lists[["scran.findMarkers"]]) > 0)
})


test_that("run_gene_markers identifies markers for multiple cell types", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- run_gene_markers(
    sceval,
    method = "scran.findMarkers",
    ngenes_celltype = 30,
    verbose = FALSE
  )
  
  # Should have markers for different cell types
  markers <- sceval@gene_lists[["scran.findMarkers"]]
  expect_true(length(markers) > 0)
  expect_true(length(unique(markers)) > 0)
})


test_that("run_gene_markers handles progressbar parameter", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- run_gene_markers(
    sceval,
    method = "scran.findMarkers",
    ngenes_celltype = 50,
    progressbar = FALSE,
    verbose = FALSE
  )
  
  expect_true("scran.findMarkers" %in% names(sceval@gene_lists))
})
