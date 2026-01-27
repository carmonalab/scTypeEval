# Test suite for parallelization capabilities
# 
# This file tests functions that support parallel processing through
# ncores and bparam parameters

# Tests for Run.HVG parallelization ----------------------------------------

test_that("Run.HVG works with ncores = 1", {
  sceval <- create_test_scTypeEval()
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    verbose = FALSE
  )
  
  result <- Run.HVG(
    sceval,
    var.method = "basic",
    ncores = 1,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true(length(result@gene.lists) > 0)
})

test_that("Run.HVG respects ncores parameter", {
  sceval <- create_test_scTypeEval()
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    verbose = FALSE
  )
  
  # Test with ncores = 1
  result_serial <- Run.HVG(
    sceval,
    var.method = "basic",
    ncores = 1,
    verbose = FALSE
  )
  
  # Test with ncores = 2 (if available)
  result_parallel <- Run.HVG(
    sceval,
    var.method = "basic",
    ncores = 2,
    verbose = FALSE
  )
  
  # Results should have same gene lists
  expect_equal(length(result_serial@gene.lists), length(result_parallel@gene.lists))
})

test_that("Run.HVG works with bparam parameter", {
  skip_if_not_installed("BiocParallel")
  
  sceval <- create_test_scTypeEval()
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    verbose = FALSE
  )
  
  # Create a serial param
  bparam <- BiocParallel::SerialParam()
  
  result <- Run.HVG(
    sceval,
    var.method = "scran",
    bparam = bparam,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true(length(result@gene.lists) > 0)
})

test_that("Run.HVG progressbar parameter doesn't break execution", {
  sceval <- create_test_scTypeEval()
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    verbose = FALSE
  )
  
  # Should not error with progressbar = TRUE
  result <- Run.HVG(
    sceval,
    var.method = "basic",
    ncores = 1,
    progressbar = TRUE,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
})

# Tests for Run.Dissimilarity parallelization (parallelized methods only) --

# Tests for Wasserstein dissimilarity parallelization ----------------------

test_that("Run.Dissimilarity WasserStein works with ncores = 1", {
  sceval <- create_processed_scTypeEval()
  
  result <- Run.Dissimilarity(
    sceval,
    method = "WasserStein",
    ncores = 1,
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true("WasserStein" %in% names(result@dissimilarity))
})

test_that("Run.Dissimilarity WasserStein respects ncores parameter", {
  sceval <- create_processed_scTypeEval()
  
  # Test with ncores = 1
  result_serial <- Run.Dissimilarity(
    sceval,
    method = "WasserStein",
    ncores = 1,
    reduction = FALSE,
    verbose = FALSE
  )
  
  # Test with ncores = 2
  result_parallel <- Run.Dissimilarity(
    sceval,
    method = "WasserStein",
    ncores = 2,
    reduction = FALSE,
    verbose = FALSE
  )
  
  # Both should have WasserStein results with identical dist matrices
  expect_true("WasserStein" %in% names(result_serial@dissimilarity))
  expect_true("WasserStein" %in% names(result_parallel@dissimilarity))
  expect_equal(
    attr(result_serial@dissimilarity[["WasserStein"]]@dissimilarity, "Size"),
    attr(result_parallel@dissimilarity[["WasserStein"]]@dissimilarity, "Size")
  )
  expect_equal(
    as.matrix(result_serial@dissimilarity[["WasserStein"]]@dissimilarity),
    as.matrix(result_parallel@dissimilarity[["WasserStein"]]@dissimilarity)
  )
})

test_that("Run.Dissimilarity WasserStein works with bparam parameter", {
  skip_if_not_installed("BiocParallel")
  
  sceval <- create_processed_scTypeEval()
  
  # Create a serial param
  bparam <- BiocParallel::SerialParam()
  
  result <- Run.Dissimilarity(
    sceval,
    method = "WasserStein",
    bparam = bparam,
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true("WasserStein" %in% names(result@dissimilarity))
})

test_that("Run.Dissimilarity WasserStein progressbar parameter works", {
  sceval <- create_processed_scTypeEval()
  
  result <- Run.Dissimilarity(
    sceval,
    method = "WasserStein",
    ncores = 1,
    progressbar = TRUE,
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true("WasserStein" %in% names(result@dissimilarity))
})

test_that("Sequential and parallel WasserStein results are consistent", {
  sceval <- create_processed_scTypeEval()
  
  # Get sequential result
  result_seq <- Run.Dissimilarity(
    sceval,
    method = "WasserStein",
    ncores = 1,
    reduction = FALSE,
    verbose = FALSE
  )
  
  # Get parallel result
  result_par <- Run.Dissimilarity(
    sceval,
    method = "WasserStein",
    ncores = 2,
    reduction = FALSE,
    verbose = FALSE
  )
  
  # Should have identical dist matrices
  expect_equal(
    attr(result_seq@dissimilarity[["WasserStein"]]@dissimilarity, "Size"),
    attr(result_par@dissimilarity[["WasserStein"]]@dissimilarity, "Size")
  )
  expect_equal(
    as.matrix(result_seq@dissimilarity[["WasserStein"]]@dissimilarity),
    as.matrix(result_par@dissimilarity[["WasserStein"]]@dissimilarity)
  )
})

test_that("Run.Dissimilarity WasserStein works with multiple cores in large dataset", {
  sceval <- create_processed_scTypeEval(small = FALSE)
  
  result <- Run.Dissimilarity(
    sceval,
    method = "WasserStein",
    ncores = 2,
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true("WasserStein" %in% names(result@dissimilarity))
})

# Tests for Run.GeneMarkers parallelization --------------------------------

test_that("Run.GeneMarkers works with ncores = 1", {
  sceval <- create_processed_scTypeEval()
  
  result <- Run.GeneMarkers(
    sceval,
    ncores = 1,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true(length(result@gene.lists) > 0)
})

test_that("Run.GeneMarkers respects ncores parameter", {
  sceval <- create_processed_scTypeEval()
  
  # Test with ncores = 1
  result_serial <- Run.GeneMarkers(
    sceval,
    ncores = 1,
    verbose = FALSE
  )
  
  # Test with ncores = 2
  result_parallel <- Run.GeneMarkers(
    sceval,
    ncores = 2,
    verbose = FALSE
  )
  
  # Both should have gene markers
  expect_true(length(result_serial@gene.lists) > 0)
  expect_true(length(result_parallel@gene.lists) > 0)
  expect_equal(result_serial@gene.lists, result_parallel@gene.lists)
})

test_that("Run.GeneMarkers works with bparam parameter", {
  skip_if_not_installed("BiocParallel")
  
  sceval <- create_processed_scTypeEval()
  
  # Create a serial param
  bparam <- BiocParallel::SerialParam()
  
  result <- Run.GeneMarkers(
    sceval,
    bparam = bparam,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true(length(result@gene.lists) > 0)
})

# Tests for RecipClassif parallelization -----------------------------------

test_that("RecipClassif works with ncores = 1", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(
    sceval,
    method = "RecipClassif:Match",
    ncores = 1,
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_s4_class(sceval, "scTypeEval")
  expect_true("RecipClassif:Match" %in% names(sceval@dissimilarity))
})

test_that("RecipClassif respects ncores parameter", {
  sceval <- create_processed_scTypeEval()
  
  # Test with ncores = 1
  result_serial <- Run.Dissimilarity(
    sceval,
    method = "RecipClassif:Match",
    ncores = 1,
    reduction = FALSE,
    verbose = FALSE
  )
  
  # Test with ncores = 2
  result_parallel <- Run.Dissimilarity(
    sceval,
    method = "RecipClassif:Match",
    ncores = 2,
    reduction = FALSE,
    verbose = FALSE
  )
  
  # Both should have identical results
  expect_true("RecipClassif:Match" %in% names(result_serial@dissimilarity))
  expect_true("RecipClassif:Match" %in% names(result_parallel@dissimilarity))
  expect_equal(
    as.matrix(result_serial@dissimilarity[["RecipClassif:Match"]]@dissimilarity),
    as.matrix(result_parallel@dissimilarity[["RecipClassif:Match"]]@dissimilarity)
  )
})

test_that("RecipClassif works with bparam parameter", {
  skip_if_not_installed("BiocParallel")
  
  sceval <- create_processed_scTypeEval()
  
  # Create a serial param
  bparam <- BiocParallel::SerialParam()
  
  result <- Run.Dissimilarity(
    sceval,
    method = "RecipClassif:Match",
    bparam = bparam,
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true("RecipClassif:Match" %in% names(result@dissimilarity))
})

# Integration tests: Sequential vs Parallel consistency -------------------

test_that("Sequential and parallel HVG results are consistent", {
  sceval <- create_test_scTypeEval()
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    verbose = FALSE
  )
  
  # Get sequential result
  result_seq <- Run.HVG(
    sceval,
    var.method = "basic",
    ncores = 1,
    verbose = FALSE
  )
  
  # Get parallel result (same seed for reproducibility)
  result_par <- Run.HVG(
    sceval,
    var.method = "basic",
    ncores = 2,
    verbose = FALSE
  )
  
  # Should have same number of gene lists
  expect_equal(result_seq@gene.lists, result_par@gene.lists)
})

test_that("Sequential and parallel Dissimilarity results are consistent", {
  sceval <- create_processed_scTypeEval()
  
  # Get sequential result
  result_seq <- Run.Dissimilarity(
    sceval,
    method = "RecipClassif:Match",
    ncores = 1,
    reduction = FALSE,
    verbose = FALSE
  )
  
  # Get parallel result
  result_par <- Run.Dissimilarity(
    sceval,
    method = "RecipClassif:Match",
    ncores = 2,
    reduction = FALSE,
    verbose = FALSE
  )
  
  # Should have identical dist matrices
  expect_equal(
    attr(result_seq@dissimilarity[["RecipClassif:Match"]]@dissimilarity, "Size"),
    attr(result_par@dissimilarity[["RecipClassif:Match"]]@dissimilarity, "Size")
  )
  expect_equal(
    as.matrix(result_seq@dissimilarity[["RecipClassif:Match"]]@dissimilarity),
    as.matrix(result_par@dissimilarity[["RecipClassif:Match"]]@dissimilarity)
  )
})

# Tests for edge cases in parallelization -----------------------------------

test_that("Parallelization works with single sample", {
  sceval <- create_test_scTypeEval()
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    min.samples = 1,  # Allow single sample
    verbose = FALSE
  )
  
  result <- Run.HVG(
    sceval,
    var.method = "basic",
    ncores = 2,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
})

test_that("Progressbar works with parallelization", {
  sceval <- create_test_scTypeEval()
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    verbose = FALSE
  )
  
  result <- Run.HVG(
    sceval,
    var.method = "basic",
    ncores = 2,
    progressbar = TRUE,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
})
