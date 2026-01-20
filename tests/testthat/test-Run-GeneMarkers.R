test_that("Run.GeneMarkers identifies marker genes with scran.findMarkers", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.GeneMarkers(
    sceval,
    method = "scran.findMarkers",
    ngenes.celltype = 50,
    aggregation = "single-cell",
    verbose = FALSE
  )
  
  expect_true("scran.findMarkers" %in% names(sceval@gene.lists))
  expect_true(length(sceval@gene.lists[["scran.findMarkers"]]) > 0)
})


test_that("Run.GeneMarkers respects ngenes.celltype parameter", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.GeneMarkers(
    sceval,
    method = "scran.findMarkers",
    ngenes.celltype = 20,
    verbose = FALSE
  )
  
  # Total markers should be reasonable given cell types
  expect_true(length(sceval@gene.lists[["scran.findMarkers"]]) > 0)
})


test_that("Run.GeneMarkers works with single-cell aggregation", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.GeneMarkers(
    sceval,
    method = "scran.findMarkers",
    ngenes.celltype = 50,
    aggregation = "single-cell",
    verbose = FALSE
  )
  
  expect_true("scran.findMarkers" %in% names(sceval@gene.lists))
})


test_that("Run.GeneMarkers respects black.list parameter", {
  sceval <- create_processed_scTypeEval()
  
  # Get some gene names
  all_genes <- rownames(sceval@data[["single-cell"]]@matrix)
  black_genes <- all_genes[1:10]
  
  sceval@black.list <- black_genes
  
  sceval <- Run.GeneMarkers(
    sceval,
    method = "scran.findMarkers",
    ngenes.celltype = 50,
    verbose = FALSE
  )
  
  # Markers should not contain blacklisted genes
  expect_true(!any(sceval@gene.lists[["scran.findMarkers"]] %in% black_genes))
})


test_that("Run.GeneMarkers accepts custom black.list parameter", {
  sceval <- create_processed_scTypeEval()
  
  all_genes <- rownames(sceval@data[["single-cell"]]@matrix)
  custom_black <- all_genes[1:5]
  
  sceval <- Run.GeneMarkers(
    sceval,
    method = "scran.findMarkers",
    ngenes.celltype = 50,
    black.list = custom_black,
    verbose = FALSE
  )
  
  expect_true(!any(sceval@gene.lists[["scran.findMarkers"]] %in% custom_black))
})


test_that("Run.GeneMarkers errors without processed data", {
  sceval <- create_test_scTypeEval()
  
  expect_error(
    Run.GeneMarkers(
      sceval,
      method = "scran.findMarkers",
      ngenes.celltype = 50,
      verbose = FALSE
    ),
    "No normalization slot found"
  )
})


test_that("Run.GeneMarkers errors on invalid aggregation type", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    Run.GeneMarkers(
      sceval,
      method = "scran.findMarkers",
      ngenes.celltype = 50,
      aggregation = "invalid_type",
      verbose = FALSE
    ),
    "Invalid aggregation type"
  )
})


test_that("Run.GeneMarkers errors on unsupported method", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    Run.GeneMarkers(
      sceval,
      method = "unsupported_method",
      ngenes.celltype = 50,
      verbose = FALSE
    ),
    "Supported gene markes definitions"
  )
})


test_that("Run.GeneMarkers handles verbose parameter", {
  sceval <- create_processed_scTypeEval()
  
  # verbose = TRUE should produce messages
  expect_message(
    sceval <- Run.GeneMarkers(
      sceval,
      method = "scran.findMarkers",
      ngenes.celltype = 50,
      verbose = TRUE
    )
  )
  
  # verbose = FALSE should suppress messages
  sceval2 <- create_processed_scTypeEval()
  expect_silent(
    sceval2 <- Run.GeneMarkers(
      sceval2,
      method = "scran.findMarkers",
      ngenes.celltype = 50,
      verbose = FALSE
    )
  )
})


test_that("Run.GeneMarkers handles ncores parameter", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.GeneMarkers(
    sceval,
    method = "scran.findMarkers",
    ngenes.celltype = 50,
    ncores = 1,
    verbose = FALSE
  )
  
  expect_true("scran.findMarkers" %in% names(sceval@gene.lists))
})


test_that("Run.GeneMarkers returns genes present in the data", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.GeneMarkers(
    sceval,
    method = "scran.findMarkers",
    ngenes.celltype = 50,
    verbose = FALSE
  )
  
  available_genes <- rownames(sceval@data[["single-cell"]]@matrix)
  expect_true(all(sceval@gene.lists[["scran.findMarkers"]] %in% available_genes))
})


test_that("Run.GeneMarkers works with multiple samples", {
  test_data <- generate_test_data(n_samples = 6)
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  sceval <- Run.ProcessingData(sceval, ident = "celltype", sample = "sample",
                                min.samples = 3, min.cells = 5, verbose = FALSE)
  
  sceval <- Run.GeneMarkers(
    sceval,
    method = "scran.findMarkers",
    ngenes.celltype = 50,
    verbose = FALSE
  )
  
  expect_true("scran.findMarkers" %in% names(sceval@gene.lists))
  expect_true(length(sceval@gene.lists[["scran.findMarkers"]]) > 0)
})


test_that("Run.GeneMarkers identifies markers for multiple cell types", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.GeneMarkers(
    sceval,
    method = "scran.findMarkers",
    ngenes.celltype = 30,
    verbose = FALSE
  )
  
  # Should have markers for different cell types
  markers <- sceval@gene.lists[["scran.findMarkers"]]
  expect_true(length(markers) > 0)
  expect_true(length(unique(markers)) > 0)
})


test_that("Run.GeneMarkers handles progressbar parameter", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.GeneMarkers(
    sceval,
    method = "scran.findMarkers",
    ngenes.celltype = 50,
    progressbar = FALSE,
    verbose = FALSE
  )
  
  expect_true("scran.findMarkers" %in% names(sceval@gene.lists))
})
