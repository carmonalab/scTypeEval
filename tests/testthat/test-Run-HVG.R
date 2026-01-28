test_that("Run.HVG identifies highly variable genes with scran method", {
  skip_if_not_installed("bluster")
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  sceval <- Run.HVG(
    sceval,
    var.method = "scran",
    ngenes = 100,
    sample = TRUE,
    aggregation = "single-cell",
    verbose = FALSE
  )
  
  expect_true("HVG" %in% names(sceval@gene.lists))
  expect_true(length(sceval@gene.lists$HVG) > 0)
  expect_true(length(sceval@gene.lists$HVG) <= 100)
})


test_that("Run.HVG identifies highly variable genes with basic method", {
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  sceval <- Run.HVG(
    sceval,
    var.method = "basic",
    ngenes = 100,
    sample = TRUE,
    aggregation = "single-cell",
    verbose = FALSE
  )
  
  expect_true("HVG" %in% names(sceval@gene.lists))
  expect_true(length(sceval@gene.lists$HVG) > 0)
})


test_that("Run.HVG basic respects ngenes parameter", {
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  sceval <- Run.HVG(
    sceval,
    var.method = "basic",
    ngenes = 50,
    sample = TRUE,
    verbose = FALSE
  )
  
  expect_true(length(sceval@gene.lists$HVG) <= 50)
})


test_that("Run.HVG basic works with sample blocking", {
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  sceval <- Run.HVG(
    sceval,
    var.method = "basic",
    ngenes = 100,
    sample = TRUE,
    verbose = FALSE
  )
  
  expect_true("HVG" %in% names(sceval@gene.lists))
})


test_that("Run.HVG basic works without sample blocking", {
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  expect_warning(
    sceval <- Run.HVG(
      sceval,
      var.method = "basic",
      ngenes = 100,
      sample = FALSE,
      verbose = FALSE
    ),
    "Not leveraging sample information"
  )
  
  expect_true("HVG" %in% names(sceval@gene.lists))
})

test_that("Run.HVG scran respects ngenes parameter", {
  skip_if_not_installed("bluster")
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  sceval <- Run.HVG(
    sceval,
    var.method = "scran",
    ngenes = 50,
    sample = TRUE,
    verbose = FALSE
  )
  
  expect_true(length(sceval@gene.lists$HVG) <= 50)
})


test_that("Run.HVG scran works with sample blocking", {
  skip_if_not_installed("bluster")
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  sceval <- Run.HVG(
    sceval,
    var.method = "scran",
    ngenes = 100,
    sample = TRUE,
    verbose = FALSE
  )
  
  expect_true("HVG" %in% names(sceval@gene.lists))
})


test_that("Run.HVG scran works without sample blocking", {
  skip_if_not_installed("bluster")
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  expect_warning(
    sceval <- Run.HVG(
      sceval,
      var.method = "scran",
      ngenes = 100,
      sample = FALSE,
      verbose = FALSE
    ),
    "Not leveraging sample information"
  )
  
  expect_true("HVG" %in% names(sceval@gene.lists))
})


test_that("Run.HVG basic works with pseudobulk aggregation", {
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  sceval <- Run.HVG(
    sceval,
    var.method = "basic",
    ngenes = 100,
    aggregation = "pseudobulk",
    verbose = FALSE
  )
  
  expect_true("HVG" %in% names(sceval@gene.lists))
})

test_that("Run.HVG scran works with pseudobulk aggregation", {
  skip_if_not_installed("bluster")
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  sceval <- Run.HVG(
    sceval,
    var.method = "scran",
    ngenes = 100,
    aggregation = "pseudobulk",
    verbose = FALSE
  )
  
  expect_true("HVG" %in% names(sceval@gene.lists))
})



test_that("Run.HVG respects black.list parameter", {
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  # Get some gene names
  all_genes <- rownames(sceval@data[["single-cell"]]@matrix)
  black_genes <- all_genes[1:10]
  
  sceval@black.list <- black_genes
  
  sceval <- Run.HVG(
    sceval,
    var.method = "basic",
    ngenes = 100,
    verbose = FALSE
  )
  
  # HVG should not contain blacklisted genes
  expect_true(!any(sceval@gene.lists$HVG %in% black_genes))
})


test_that("Run.HVG accepts custom black.list parameter", {
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  all_genes <- rownames(sceval@data[["single-cell"]]@matrix)
  custom_black <- all_genes[1:5]
  
  sceval <- Run.HVG(
    sceval,
    var.method = "basic",
    ngenes = 100,
    black.list = custom_black,
    verbose = FALSE
  )
  
  expect_true(!any(sceval@gene.lists$HVG %in% custom_black))
})


test_that("Run.HVG errors without processed data", {
  sceval <- create_test_scTypeEval()
  
  expect_error(
    Run.HVG(
      sceval,
      var.method = "scran",
      ngenes = 100,
      verbose = FALSE
    ),
    "No normalization slot found"
  )
})


test_that("Run.HVG errors on invalid aggregation type", {
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  expect_error(
    Run.HVG(
      sceval,
      var.method = "scran",
      ngenes = 100,
      aggregation = "invalid_type",
      verbose = FALSE
    ),
    "Invalid aggregation type"
  )
})


test_that("Run.HVG errors on unsupported var.method", {
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  expect_error(
    Run.HVG(
      sceval,
      var.method = "unsupported_method",
      ngenes = 100,
      verbose = FALSE
    ),
    "not supported for getting variable genes"
  )
})


test_that("Run.HVG handles verbose parameter", {
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  # verbose = TRUE should produce messages
  expect_message(
    sceval <- Run.HVG(
      sceval,
      var.method = "basic",
      ngenes = 100,
      verbose = TRUE
    ),
    "Computing HVG"
  )
  
  # verbose = FALSE should suppress messages
  sceval2 <- create_processed_scTypeEval(HVG = FALSE)
  expect_silent(
    sceval2 <- Run.HVG(
      sceval2,
      var.method = "basic",
      ngenes = 100,
      verbose = FALSE
    )
  )
})


test_that("Run.HVG handles ncores parameter", {
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  sceval <- Run.HVG(
    sceval,
    var.method = "basic",
    ngenes = 100,
    ncores = 1,
    verbose = FALSE
  )
  
  expect_true("HVG" %in% names(sceval@gene.lists))
})


test_that("Run.HVG returns genes present in the data", {
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  sceval <- Run.HVG(
    sceval,
    var.method = "basic",
    ngenes = 100,
    verbose = FALSE
  )
  
  available_genes <- rownames(sceval@data[["single-cell"]]@matrix)
  expect_true(all(sceval@gene.lists$HVG %in% available_genes))
})

