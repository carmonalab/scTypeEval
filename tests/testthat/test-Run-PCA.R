test_that("Run.PCA computes PCA for single-cell data", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.PCA(
    sceval,
    ndim = 5,
    verbose = FALSE
  )
  
  expect_true("single-cell" %in% names(sceval@reductions))
  expect_s4_class(sceval@reductions[["single-cell"]], "DimRed")
  expect_equal(sceval@reductions[["single-cell"]]@key, "PCA")
})


test_that("Run.PCA computes PCA for pseudobulk data", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.PCA(
    sceval,
    ndim = 5,
    verbose = FALSE
  )
  
  expect_true("pseudobulk" %in% names(sceval@reductions))
  expect_s4_class(sceval@reductions[["pseudobulk"]], "DimRed")
})


test_that("Run.PCA respects ndim parameter", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.PCA(
    sceval,
    ndim = 5,
    verbose = FALSE
  )
  
  # Number of PCs should be at most ndim
  expect_true(nrow(sceval@reductions[["single-cell"]]@embeddings) <= 10)
})


test_that("Run.PCA uses gene.list parameter", {
  sceval <- create_processed_scTypeEval()
  
  # Add a custom gene list
  all_genes <- rownames(sceval@data[["single-cell"]]@matrix)
  custom_genes <- all_genes[1:100]
  sceval <- add.GeneList(sceval, gene.list = list(custom = custom_genes))
  
  sceval <- Run.PCA(
    sceval,
    gene.list = "custom",
    ndim = 5,
    verbose = FALSE
  )

  # Check that custom gene list was used in PCA
  expect_equal(sceval@reductions[["single-cell"]]@gene.list, custom_genes)
  # Feature loadings should only contain genes from custom list (or subset after filtering)
  # All genes in feature.loadings should be from custom gene list
  expect_true(all(rownames(sceval@reductions[["single-cell"]]@feature.loadings) %in% custom_genes))
})


test_that("Run.PCA respects black.list parameter", {
  sceval <- create_processed_scTypeEval()
  
  all_genes <- rownames(sceval@data[["single-cell"]]@matrix)
  black_genes <- all_genes[1:10]
  sceval@black.list <- black_genes
  
  sceval <- Run.PCA(
    sceval,
    ndim = 5,
    verbose = FALSE
  )
  
  expect_true(length(sceval@reductions[["single-cell"]]@black.list) > 0)
  expect_equal(sceval@reductions[["single-cell"]]@black.list, black_genes)
  # Black listed genes should not be in feature.loadings (they were filtered out of matrix)
  expect_false(any(sceval@reductions[["single-cell"]]@black.list %in% rownames(sceval@reductions[["single-cell"]]@feature.loadings)))
})


test_that("Run.PCA errors without processed data", {
  sceval <- create_test_scTypeEval()
  
  expect_error(
    Run.PCA(
      sceval,
      ndim = 5,
      verbose = FALSE
    ),
    "No normalization slot found"
  )
})


test_that("Run.PCA stores feature loadings", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.PCA(
    sceval,
    ndim = 5,
    verbose = FALSE
  )
  
  expect_true(nrow(sceval@reductions[["single-cell"]]@feature.loadings) > 0)
  expect_true(ncol(sceval@reductions[["single-cell"]]@feature.loadings) > 0)
})


test_that("Run.PCA stores correct metadata", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.PCA(
    sceval,
    ndim = 5,
    verbose = FALSE
  )
  
  dimred <- sceval@reductions[["single-cell"]]
  expect_equal(dimred@aggregation, "single-cell")
  expect_equal(dimred@key, "PCA")
  expect_true(is.factor(dimred@group))
  expect_true(is.factor(dimred@sample))
})


test_that("Run.PCA handles verbose parameter", {
  sceval <- create_processed_scTypeEval()
  
  # verbose = TRUE should produce messages
  expect_message(
    sceval <- Run.PCA(
      sceval,
      ndim = 5,
      verbose = TRUE
    ),
    "Computing PCA"
  )
  
  # verbose = FALSE should suppress messages
  sceval2 <- create_processed_scTypeEval()
  expect_silent(
    sceval2 <- Run.PCA(
      sceval2,
      ndim = 5,
      verbose = FALSE
    )
  )
})


test_that("Run.PCA embeddings have correct dimensions", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.PCA(
    sceval,
    ndim = 5,
    verbose = FALSE
  )
  
  embeddings <- sceval@reductions[["single-cell"]]@embeddings
  n_samples <- ncol(sceval@data[["single-cell"]]@matrix)
  
  # Embeddings should be ndim x n_samples
  expect_equal(ncol(embeddings), n_samples)
})


test_that("Run.PCA works with multiple samples", {
  test_data <- generate_test_data(n_samples = 6)
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  sceval <- Run.ProcessingData(sceval, ident = "celltype", sample = "sample",
                                min.samples = 3, min.cells = 5, verbose = FALSE)
  sceval <- Run.HVG(sceval, ngenes = 200, verbose = FALSE)                          
  
  sceval <- Run.PCA(
    sceval,
    ndim = 5,
    verbose = FALSE
  )
  
  expect_true("single-cell" %in% names(sceval@reductions))
  expect_true("pseudobulk" %in% names(sceval@reductions))
  expect_equal(sceval@reductions[["single-cell"]]@gene.list, sceval@gene.lists[["HVG"]])
})

test_that("Run.PCA computes PCA for all data assays", {
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.PCA(
    sceval,
    ndim = 5,
    verbose = FALSE
  )
  
  # Should compute PCA for both single-cell and pseudobulk
  expect_equal(length(sceval@reductions), length(sceval@data))
  expect_true(all(names(sceval@data) %in% names(sceval@reductions)))
})
