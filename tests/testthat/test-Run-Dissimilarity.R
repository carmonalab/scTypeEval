test_that("run_dissimilarity computes Pseudobulk:Euclidean dissimilarity", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  sceval <- run_dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_true("Pseudobulk:Euclidean" %in% names(sceval@dissimilarity))
  expect_s4_class(sceval@dissimilarity[["Pseudobulk:Euclidean"]], "dissimilarity_assay")
  expect_s3_class(sceval@dissimilarity[["Pseudobulk:Euclidean"]]@dissimilarity, "dist")
})


test_that("run_dissimilarity computes Pseudobulk:Cosine dissimilarity", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  sceval <- run_dissimilarity(
    sceval,
    method = "Pseudobulk:Cosine",
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_true("Pseudobulk:Cosine" %in% names(sceval@dissimilarity))
  expect_equal(sceval@dissimilarity[["Pseudobulk:Cosine"]]@method, "Pseudobulk:Cosine")
})


test_that("run_dissimilarity computes Pseudobulk:Pearson dissimilarity", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  sceval <- run_dissimilarity(
    sceval,
    method = "Pseudobulk:Pearson",
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_true("Pseudobulk:Pearson" %in% names(sceval@dissimilarity))
  expect_equal(sceval@dissimilarity[["Pseudobulk:Pearson"]]@method, "Pseudobulk:Pearson")
})


test_that("run_dissimilarity computes WasserStein dissimilarity", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  sceval <- run_dissimilarity(
    sceval,
    method = "WasserStein",
    reduction = TRUE,
    verbose = FALSE
  )
  
  expect_true("WasserStein" %in% names(sceval@dissimilarity))
  expect_equal(sceval@dissimilarity[["WasserStein"]]@method, "WasserStein")
})


test_that("run_dissimilarity computes recip_classif:Match dissimilarity", {
  skip_if_not_installed("SingleR")
  
  sceval <- create_processed_scTypeEval()
  
  sceval <- run_dissimilarity(
    sceval,
    method = "recip_classif:Match",
    reciprocal_classifier = "SingleR",
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_true("recip_classif:Match" %in% names(sceval@dissimilarity))
  expect_equal(sceval@dissimilarity[["recip_classif:Match"]]@method, "recip_classif:Match")
})


test_that("run_dissimilarity computes recip_classif:Score dissimilarity", {
  skip_if_not_installed("SingleR")
  
  sceval <- create_processed_scTypeEval()
  
  sceval <- run_dissimilarity(
    sceval,
    method = "recip_classif:Score",
    reciprocal_classifier = "SingleR",
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_true("recip_classif:Score" %in% names(sceval@dissimilarity))
  expect_equal(sceval@dissimilarity[["recip_classif:Score"]]@method, "recip_classif:Score")
})


test_that("run_dissimilarity works with reduction = TRUE", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  sceval <- run_dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = TRUE,
    verbose = FALSE
  )
  
  expect_true("Pseudobulk:Euclidean" %in% names(sceval@dissimilarity))
})


test_that("run_dissimilarity works with reduction = FALSE", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- run_dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_true("Pseudobulk:Euclidean" %in% names(sceval@dissimilarity))
})


test_that("run_dissimilarity uses gene_list parameter", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- run_dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = FALSE,
    gene_list = "HVG",
    verbose = FALSE
  )
  
  expect_equal(sceval@dissimilarity[["Pseudobulk:Euclidean"]]@gene_list, sceval@gene_lists[["HVG"]])
})


test_that("run_dissimilarity respects black_list parameter", {
  sceval <- create_processed_scTypeEval()
  
  all_genes <- rownames(sceval@data[["pseudobulk"]]@matrix)
  black_genes <- all_genes[1:10]
  
  sceval <- run_dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = FALSE,
    black_list = black_genes,
    verbose = FALSE
  )
  
  expect_true(length(sceval@dissimilarity[["Pseudobulk:Euclidean"]]@black_list) > 0)
})


test_that("run_dissimilarity errors without processed data", {
  sceval <- create_test_scTypeEval()
  
  expect_error(
    run_dissimilarity(
      sceval,
      method = "Pseudobulk:Euclidean",
      reduction = FALSE,
      verbose = FALSE
    ),
    "No processed data slot found"
  )
})


test_that("run_dissimilarity errors without PCA when reduction = TRUE", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    run_dissimilarity(
      sceval,
      method = "Pseudobulk:Euclidean",
      reduction = TRUE,
      verbose = FALSE
    ),
    "No dimensional reduction slot found"
  )
})


test_that("run_dissimilarity errors on unsupported method", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    run_dissimilarity(
      sceval,
      method = "UnsupportedMethod",
      reduction = FALSE,
      verbose = FALSE
    ),
    "Not supported dissimilarity method"
  )
})


test_that("run_dissimilarity handles verbose parameter", {
  sceval <- create_processed_scTypeEval()
  
  # verbose = TRUE should produce messages
  expect_message(
    sceval <- run_dissimilarity(
      sceval,
      method = "Pseudobulk:Euclidean",
      reduction = FALSE,
      verbose = TRUE
    )
  )
  
  # verbose = FALSE should suppress messages
  sceval2 <- create_processed_scTypeEval()
  expect_silent(
    sceval2 <- run_dissimilarity(
      sceval2,
      method = "Pseudobulk:Euclidean",
      reduction = FALSE,
      verbose = FALSE
    )
  )
})


test_that("run_dissimilarity handles ncores parameter", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- run_dissimilarity(
    sceval,
    method = "WasserStein",
    reduction = FALSE,
    ncores = 1,
    verbose = FALSE
  )
  
  expect_true("WasserStein" %in% names(sceval@dissimilarity))
})


test_that("run_dissimilarity stores correct metadata", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- run_dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = FALSE,
    verbose = FALSE
  )
  
  diss_assay <- sceval@dissimilarity[["Pseudobulk:Euclidean"]]
  expect_equal(diss_assay@method, "Pseudobulk:Euclidean")
  expect_equal(diss_assay@aggregation, "pseudobulk")
  expect_true(is.list(diss_assay@ident))
  expect_true(is.factor(diss_assay@sample))
})

test_that("run_dissimilarity handles progressbar parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  sceval <- run_dissimilarity(
    sceval,
    method = "WasserStein",
    reduction = TRUE,
    progressbar = FALSE,
    verbose = FALSE
  )
  
  expect_true("WasserStein" %in% names(sceval@dissimilarity))
})


test_that("run_dissimilarity switches reduction to FALSE for recip_classif methods", {
  skip_if_not_installed("SingleR")
  
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  expect_warning(
    sceval <- run_dissimilarity(
      sceval,
      method = "recip_classif:Match",
      reciprocal_classifier = "SingleR",
      reduction = TRUE,
      verbose = FALSE
    ),
    "No dimensional reduction dissimilarity computation supported"
  )
})


test_that("run_dissimilarity uses multiple cores for WasserStein method", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  # Test with ncores = 2
  sceval <- run_dissimilarity(
    sceval,
    method = "WasserStein",
    reduction = TRUE,
    ncores = 2,
    verbose = FALSE
  )
  
  expect_true("WasserStein" %in% names(sceval@dissimilarity))
  expect_s4_class(sceval@dissimilarity[["WasserStein"]], "dissimilarity_assay")
  expect_equal(sceval@dissimilarity[["WasserStein"]]@method, "WasserStein")
})


test_that("run_dissimilarity uses multiple cores for recip_classif methods", {
  skip_if_not_installed("SingleR")
  
  sceval <- create_processed_scTypeEval()
  
  # Test recip_classif:Match with ncores = 2
  sceval <- run_dissimilarity(
    sceval,
    method = "recip_classif:Match",
    reciprocal_classifier = "SingleR",
    reduction = FALSE,
    ncores = 2,
    verbose = FALSE
  )
  
  expect_true("recip_classif:Match" %in% names(sceval@dissimilarity))
  expect_equal(sceval@dissimilarity[["recip_classif:Match"]]@method, "recip_classif:Match")
})


test_that("run_dissimilarity ncores parameter works for both parallelizable methods", {
  skip_if_not_installed("SingleR")
  
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  # Test WasserStein
  sceval <- run_dissimilarity(
    sceval,
    method = "WasserStein",
    reduction = TRUE,
    ncores = 2,
    verbose = FALSE
  )
  
  # Test recip_classif:Score
  sceval <- run_dissimilarity(
    sceval,
    method = "recip_classif:Score",
    reciprocal_classifier = "SingleR",
    reduction = FALSE,
    ncores = 2,
    verbose = FALSE
  )
  
  expect_true("WasserStein" %in% names(sceval@dissimilarity))
  expect_true("recip_classif:Score" %in% names(sceval@dissimilarity))
})


test_that("run_dissimilarity ncores parameter accepted but not used for Pseudobulk methods", {
  sceval <- create_processed_scTypeEval()
  
  # Pseudobulk methods don't use parallelization, but ncores parameter should be accepted
  sceval <- run_dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = FALSE,
    ncores = 2,
    verbose = FALSE
  )
  
  expect_true("Pseudobulk:Euclidean" %in% names(sceval@dissimilarity))
  
  sceval <- run_dissimilarity(
    sceval,
    method = "Pseudobulk:Cosine",
    reduction = FALSE,
    ncores = 2,
    verbose = FALSE
  )
  
  expect_true("Pseudobulk:Cosine" %in% names(sceval@dissimilarity))
})
