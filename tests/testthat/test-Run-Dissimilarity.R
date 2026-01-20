test_that("Run.Dissimilarity computes Pseudobulk:Euclidean dissimilarity", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  sceval <- Run.Dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_true("Pseudobulk:Euclidean" %in% names(sceval@dissimilarity))
  expect_s4_class(sceval@dissimilarity[["Pseudobulk:Euclidean"]], "DissimilarityAssay")
  expect_s3_class(sceval@dissimilarity[["Pseudobulk:Euclidean"]]@dissimilarity, "dist")
})


test_that("Run.Dissimilarity computes Pseudobulk:Cosine dissimilarity", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  sceval <- Run.Dissimilarity(
    sceval,
    method = "Pseudobulk:Cosine",
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_true("Pseudobulk:Cosine" %in% names(sceval@dissimilarity))
  expect_equal(sceval@dissimilarity[["Pseudobulk:Cosine"]]@method, "Pseudobulk:Cosine")
})


test_that("Run.Dissimilarity computes Pseudobulk:Pearson dissimilarity", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  sceval <- Run.Dissimilarity(
    sceval,
    method = "Pseudobulk:Pearson",
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_true("Pseudobulk:Pearson" %in% names(sceval@dissimilarity))
  expect_equal(sceval@dissimilarity[["Pseudobulk:Pearson"]]@method, "Pseudobulk:Pearson")
})


test_that("Run.Dissimilarity computes WasserStein dissimilarity", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  sceval <- Run.Dissimilarity(
    sceval,
    method = "WasserStein",
    reduction = TRUE,
    verbose = FALSE
  )
  
  expect_true("WasserStein" %in% names(sceval@dissimilarity))
  expect_equal(sceval@dissimilarity[["WasserStein"]]@method, "WasserStein")
})


test_that("Run.Dissimilarity computes RecipClassif:Match dissimilarity", {
  skip_if_not_installed("SingleR")
  
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.Dissimilarity(
    sceval,
    method = "RecipClassif:Match",
    ReciprocalClassifier = "SingleR",
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_true("RecipClassif:Match" %in% names(sceval@dissimilarity))
  expect_equal(sceval@dissimilarity[["RecipClassif:Match"]]@method, "RecipClassif:Match")
})


test_that("Run.Dissimilarity computes RecipClassif:Score dissimilarity", {
  skip_if_not_installed("SingleR")
  
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.Dissimilarity(
    sceval,
    method = "RecipClassif:Score",
    ReciprocalClassifier = "SingleR",
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_true("RecipClassif:Score" %in% names(sceval@dissimilarity))
  expect_equal(sceval@dissimilarity[["RecipClassif:Score"]]@method, "RecipClassif:Score")
})


test_that("Run.Dissimilarity works with reduction = TRUE", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  sceval <- Run.Dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = TRUE,
    verbose = FALSE
  )
  
  expect_true("Pseudobulk:Euclidean" %in% names(sceval@dissimilarity))
})


test_that("Run.Dissimilarity works with reduction = FALSE", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.Dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_true("Pseudobulk:Euclidean" %in% names(sceval@dissimilarity))
})


test_that("Run.Dissimilarity uses gene.list parameter", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.Dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = FALSE,
    gene.list = "HVG",
    verbose = FALSE
  )
  
  expect_equal(sceval@dissimilarity[["Pseudobulk:Euclidean"]]@gene.list, sceval@gene.lists[["HVG"]])
})


test_that("Run.Dissimilarity respects black.list parameter", {
  sceval <- create_processed_scTypeEval()
  
  all_genes <- rownames(sceval@data[["pseudobulk"]]@matrix)
  black_genes <- all_genes[1:10]
  
  sceval <- Run.Dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = FALSE,
    black.list = black_genes,
    verbose = FALSE
  )
  
  expect_true(length(sceval@dissimilarity[["Pseudobulk:Euclidean"]]@black.list) > 0)
})


test_that("Run.Dissimilarity errors without processed data", {
  sceval <- create_test_scTypeEval()
  
  expect_error(
    Run.Dissimilarity(
      sceval,
      method = "Pseudobulk:Euclidean",
      reduction = FALSE,
      verbose = FALSE
    ),
    "No processed data slot found"
  )
})


test_that("Run.Dissimilarity errors without PCA when reduction = TRUE", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    Run.Dissimilarity(
      sceval,
      method = "Pseudobulk:Euclidean",
      reduction = TRUE,
      verbose = FALSE
    ),
    "No dimensional reduction slot found"
  )
})


test_that("Run.Dissimilarity errors on unsupported method", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    Run.Dissimilarity(
      sceval,
      method = "UnsupportedMethod",
      reduction = FALSE,
      verbose = FALSE
    ),
    "Not supported dissimilarity method"
  )
})


test_that("Run.Dissimilarity handles verbose parameter", {
  sceval <- create_processed_scTypeEval()
  
  # verbose = TRUE should produce messages
  expect_message(
    sceval <- Run.Dissimilarity(
      sceval,
      method = "Pseudobulk:Euclidean",
      reduction = FALSE,
      verbose = TRUE
    )
  )
  
  # verbose = FALSE should suppress messages
  sceval2 <- create_processed_scTypeEval()
  expect_silent(
    sceval2 <- Run.Dissimilarity(
      sceval2,
      method = "Pseudobulk:Euclidean",
      reduction = FALSE,
      verbose = FALSE
    )
  )
})


test_that("Run.Dissimilarity handles ncores parameter", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.Dissimilarity(
    sceval,
    method = "WasserStein",
    reduction = FALSE,
    ncores = 1,
    verbose = FALSE
  )
  
  expect_true("WasserStein" %in% names(sceval@dissimilarity))
})


test_that("Run.Dissimilarity stores correct metadata", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.Dissimilarity(
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


test_that("Run.Dissimilarity works with multiple samples", {
  test_data <- generate_test_data(n_samples = 6)
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  sceval <- Run.ProcessingData(sceval, ident = "celltype", sample = "sample",
                                min.samples = 3, min.cells = 5, verbose = FALSE)
  sceval <- Run.HVG(sceval, ngenes = 100, verbose = FALSE)                              
  
  sceval <- Run.Dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = FALSE,
    verbose = FALSE
  )
  
  expect_true("Pseudobulk:Euclidean" %in% names(sceval@dissimilarity))
})


test_that("Run.Dissimilarity handles progressbar parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  sceval <- Run.Dissimilarity(
    sceval,
    method = "WasserStein",
    reduction = TRUE,
    progressbar = FALSE,
    verbose = FALSE
  )
  
  expect_true("WasserStein" %in% names(sceval@dissimilarity))
})


test_that("Run.Dissimilarity switches reduction to FALSE for RecipClassif methods", {
  skip_if_not_installed("SingleR")
  
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  expect_warning(
    sceval <- Run.Dissimilarity(
      sceval,
      method = "RecipClassif:Match",
      ReciprocalClassifier = "SingleR",
      reduction = TRUE,
      verbose = FALSE
    ),
    "No dimensional reduction dissimilarity computation supported"
  )
})


test_that("Run.Dissimilarity uses multiple cores for WasserStein method", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  # Test with ncores = 2
  sceval <- Run.Dissimilarity(
    sceval,
    method = "WasserStein",
    reduction = TRUE,
    ncores = 2,
    verbose = FALSE
  )
  
  expect_true("WasserStein" %in% names(sceval@dissimilarity))
  expect_s4_class(sceval@dissimilarity[["WasserStein"]], "DissimilarityAssay")
  expect_equal(sceval@dissimilarity[["WasserStein"]]@method, "WasserStein")
})


test_that("Run.Dissimilarity uses multiple cores for RecipClassif methods", {
  skip_if_not_installed("SingleR")
  
  sceval <- create_processed_scTypeEval()
  
  # Test RecipClassif:Match with ncores = 2
  sceval <- Run.Dissimilarity(
    sceval,
    method = "RecipClassif:Match",
    ReciprocalClassifier = "SingleR",
    reduction = FALSE,
    ncores = 2,
    verbose = FALSE
  )
  
  expect_true("RecipClassif:Match" %in% names(sceval@dissimilarity))
  expect_equal(sceval@dissimilarity[["RecipClassif:Match"]]@method, "RecipClassif:Match")
})


test_that("Run.Dissimilarity ncores parameter works for both parallelizable methods", {
  skip_if_not_installed("SingleR")
  
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  # Test WasserStein
  sceval <- Run.Dissimilarity(
    sceval,
    method = "WasserStein",
    reduction = TRUE,
    ncores = 2,
    verbose = FALSE
  )
  
  # Test RecipClassif:Score
  sceval <- Run.Dissimilarity(
    sceval,
    method = "RecipClassif:Score",
    ReciprocalClassifier = "SingleR",
    reduction = FALSE,
    ncores = 2,
    verbose = FALSE
  )
  
  expect_true("WasserStein" %in% names(sceval@dissimilarity))
  expect_true("RecipClassif:Score" %in% names(sceval@dissimilarity))
})


test_that("Run.Dissimilarity ncores parameter accepted but not used for Pseudobulk methods", {
  sceval <- create_processed_scTypeEval()
  
  # Pseudobulk methods don't use parallelization, but ncores parameter should be accepted
  sceval <- Run.Dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = FALSE,
    ncores = 2,
    verbose = FALSE
  )
  
  expect_true("Pseudobulk:Euclidean" %in% names(sceval@dissimilarity))
  
  sceval <- Run.Dissimilarity(
    sceval,
    method = "Pseudobulk:Cosine",
    reduction = FALSE,
    ncores = 2,
    verbose = FALSE
  )
  
  expect_true("Pseudobulk:Cosine" %in% names(sceval@dissimilarity))
})
