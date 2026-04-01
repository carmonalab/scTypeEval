# Test suite for parallelization capabilities
# 
# This file tests functions that support parallel processing through
# ncores and bparam parameters

# Tests for run_hvg parallelization ----------------------------------------

test_that("run_hvg works with ncores = 1", {
  sceval <- create_test_scTypeEval()
  sceval <- run_processing_data(
    sceval,
    ident = "celltype",
    sample = "sample",
    verbose = FALSE
  )
  
  result <- run_hvg(
    sceval,
    var_method = "basic",
    ncores = 1,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true(length(result@gene_lists) > 0)
})

test_that("run_hvg basic respects ncores parameter", {
  sceval <- create_test_scTypeEval()
  sceval <- run_processing_data(
    sceval,
    ident = "celltype",
    sample = "sample",
    verbose = FALSE
  )
  
  # Test with ncores = 1
  result_serial <- run_hvg(
    sceval,
    var_method = "basic",
    ncores = 1,
    verbose = FALSE
  )
  
  # Test with ncores = 2 (if available)
  result_parallel <- run_hvg(
    sceval,
    var_method = "basic",
    ncores = 2,
    verbose = FALSE
  )
  
  # Results should have same gene lists
  expect_equal(length(result_serial@gene_lists), length(result_parallel@gene_lists))
})

test_that("run_hvg scran respects ncores parameter", {
  skip_if_not_installed("bluster")
  sceval <- create_test_scTypeEval()
  sceval <- run_processing_data(
    sceval,
    ident = "celltype",
    sample = "sample",
    verbose = FALSE
  )
  
  # Test with ncores = 1
  result_serial <- run_hvg(
    sceval,
    var_method = "scran",
    ncores = 1,
    verbose = FALSE
  )
  
  # Test with ncores = 2 (if available)
  result_parallel <- run_hvg(
    sceval,
    var_method = "scran",
    ncores = 2,
    verbose = FALSE
  )
  
  # Results should have same gene lists
  expect_equal(length(result_serial@gene_lists), length(result_parallel@gene_lists))
})

test_that("run_hvg works with bparam parameter", {
  skip_if_not_installed("BiocParallel")
  
  sceval <- create_test_scTypeEval()
  sceval <- run_processing_data(
    sceval,
    ident = "celltype",
    sample = "sample",
    verbose = FALSE
  )
  
  # Create a serial param
  bparam <- BiocParallel::SerialParam()
  
  result <- run_hvg(
    sceval,
    var_method = "basic",
    bparam = bparam,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true(length(result@gene_lists) > 0)
})

test_that("run_hvg progressbar parameter doesn't break execution", {
  sceval <- create_test_scTypeEval()
  sceval <- run_processing_data(
    sceval,
    ident = "celltype",
    sample = "sample",
    verbose = FALSE
  )
  
  # Should not error with progressbar = TRUE
  result <- run_hvg(
    sceval,
    var_method = "basic",
    ncores = 1,
    progressbar = TRUE,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
})

# Tests for run_dissimilarity parallelization (parallelized methods only) --

# Tests for Wasserstein dissimilarity parallelization ----------------------

test_that("run_dissimilarity WasserStein works with ncores = 1", {
  sceval <- create_processed_scTypeEval()
  
  result <- run_dissimilarity(
    sceval,
    method = "WasserStein",
    ncores = 1,
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true("WasserStein" %in% names(result@dissimilarity))
})

test_that("run_dissimilarity WasserStein respects ncores parameter", {
  sceval <- create_processed_scTypeEval()
  
  # Test with ncores = 1
  result_serial <- run_dissimilarity(
    sceval,
    method = "WasserStein",
    ncores = 1,
    reduction = FALSE,
    verbose = FALSE
  )
  
  # Test with ncores = 2
  result_parallel <- run_dissimilarity(
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

test_that("run_dissimilarity WasserStein works with bparam parameter", {
  skip_if_not_installed("BiocParallel")
  
  sceval <- create_processed_scTypeEval()
  
  # Create a serial param
  bparam <- BiocParallel::SerialParam()
  
  result <- run_dissimilarity(
    sceval,
    method = "WasserStein",
    bparam = bparam,
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true("WasserStein" %in% names(result@dissimilarity))
})

test_that("run_dissimilarity WasserStein progressbar parameter works", {
  sceval <- create_processed_scTypeEval()
  
  result <- run_dissimilarity(
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
  result_seq <- run_dissimilarity(
    sceval,
    method = "WasserStein",
    ncores = 1,
    reduction = FALSE,
    verbose = FALSE
  )
  
  # Get parallel result
  result_par <- run_dissimilarity(
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

test_that("run_dissimilarity WasserStein works with multiple cores in large dataset", {
  sceval <- create_processed_scTypeEval(small = FALSE)
  
  result <- run_dissimilarity(
    sceval,
    method = "WasserStein",
    ncores = 2,
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true("WasserStein" %in% names(result@dissimilarity))
})

# Tests for run_gene_markers parallelization --------------------------------

test_that("run_gene_markers works with ncores = 1", {
  sceval <- create_processed_scTypeEval()
  
  result <- run_gene_markers(
    sceval,
    ncores = 1,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true(length(result@gene_lists) > 0)
})

test_that("run_gene_markers respects ncores parameter", {
  sceval <- create_processed_scTypeEval()
  
  # Test with ncores = 1
  result_serial <- run_gene_markers(
    sceval,
    ncores = 1,
    verbose = FALSE
  )
  
  # Test with ncores = 2
  result_parallel <- run_gene_markers(
    sceval,
    ncores = 2,
    verbose = FALSE
  )
  
  # Both should have gene markers
  expect_true(length(result_serial@gene_lists) > 0)
  expect_true(length(result_parallel@gene_lists) > 0)
  expect_equal(result_serial@gene_lists, result_parallel@gene_lists)
})

test_that("run_gene_markers works with bparam parameter", {
  skip_if_not_installed("BiocParallel")
  
  sceval <- create_processed_scTypeEval()
  
  # Create a serial param
  bparam <- BiocParallel::SerialParam()
  
  result <- run_gene_markers(
    sceval,
    bparam = bparam,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true(length(result@gene_lists) > 0)
})

# Tests for recip_classif parallelization -----------------------------------

test_that("recip_classif works with ncores = 1", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(
    sceval,
    method = "recip_classif:Match",
    ncores = 1,
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_s4_class(sceval, "scTypeEval")
  expect_true("recip_classif:Match" %in% names(sceval@dissimilarity))
})

test_that("recip_classif respects ncores parameter", {
  sceval <- create_processed_scTypeEval()
  
  # Test with ncores = 1
  result_serial <- run_dissimilarity(
    sceval,
    method = "recip_classif:Match",
    ncores = 1,
    reduction = FALSE,
    verbose = FALSE
  )
  
  # Test with ncores = 2
  result_parallel <- run_dissimilarity(
    sceval,
    method = "recip_classif:Match",
    ncores = 2,
    reduction = FALSE,
    verbose = FALSE
  )
  
  # Both should have identical results
  expect_true("recip_classif:Match" %in% names(result_serial@dissimilarity))
  expect_true("recip_classif:Match" %in% names(result_parallel@dissimilarity))
  expect_equal(
    as.matrix(result_serial@dissimilarity[["recip_classif:Match"]]@dissimilarity),
    as.matrix(result_parallel@dissimilarity[["recip_classif:Match"]]@dissimilarity)
  )
})

test_that("recip_classif works with bparam parameter", {
  skip_if_not_installed("BiocParallel")
  
  sceval <- create_processed_scTypeEval()
  
  # Create a serial param
  bparam <- BiocParallel::SerialParam()
  
  result <- run_dissimilarity(
    sceval,
    method = "recip_classif:Match",
    bparam = bparam,
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
  expect_true("recip_classif:Match" %in% names(result@dissimilarity))
})

# Integration tests: Sequential vs Parallel consistency -------------------

test_that("Sequential and parallel HVG results are consistent", {
  sceval <- create_test_scTypeEval()
  sceval <- run_processing_data(
    sceval,
    ident = "celltype",
    sample = "sample",
    verbose = FALSE
  )
  
  # Get sequential result
  result_seq <- run_hvg(
    sceval,
    var_method = "basic",
    ncores = 1,
    verbose = FALSE
  )
  
  # Get parallel result (same seed for reproducibility)
  result_par <- run_hvg(
    sceval,
    var_method = "basic",
    ncores = 2,
    verbose = FALSE
  )
  
  # Should have same number of gene lists
  expect_equal(result_seq@gene_lists, result_par@gene_lists)
})

test_that("Sequential and parallel Dissimilarity results are consistent", {
  sceval <- create_processed_scTypeEval()
  
  # Get sequential result
  result_seq <- run_dissimilarity(
    sceval,
    method = "recip_classif:Match",
    ncores = 1,
    reduction = FALSE,
    verbose = FALSE
  )
  
  # Get parallel result
  result_par <- run_dissimilarity(
    sceval,
    method = "recip_classif:Match",
    ncores = 2,
    reduction = FALSE,
    verbose = FALSE
  )
  
  # Should have identical dist matrices
  expect_equal(
    attr(result_seq@dissimilarity[["recip_classif:Match"]]@dissimilarity, "Size"),
    attr(result_par@dissimilarity[["recip_classif:Match"]]@dissimilarity, "Size")
  )
  expect_equal(
    as.matrix(result_seq@dissimilarity[["recip_classif:Match"]]@dissimilarity),
    as.matrix(result_par@dissimilarity[["recip_classif:Match"]]@dissimilarity)
  )
})

# Tests for edge cases in parallelization -----------------------------------

test_that("Parallelization works with single sample", {
  sceval <- create_test_scTypeEval()
  sceval <- run_processing_data(
    sceval,
    ident = "celltype",
    sample = "sample",
    min_samples = 1,  # Allow single sample
    verbose = FALSE
  )
  
  result <- run_hvg(
    sceval,
    var_method = "basic",
    ncores = 2,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
})

test_that("Progressbar works with parallelization", {
  sceval <- create_test_scTypeEval()
  sceval <- run_processing_data(
    sceval,
    ident = "celltype",
    sample = "sample",
    verbose = FALSE
  )
  
  result <- run_hvg(
    sceval,
    var_method = "basic",
    ncores = 2,
    progressbar = TRUE,
    verbose = FALSE
  )
  
  expect_s4_class(result, "scTypeEval")
})
