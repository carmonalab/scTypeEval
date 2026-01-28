test_that("scTypeEval S4 class is defined correctly", {
  expect_s4_class(methods::new("scTypeEval"), "scTypeEval")
})


test_that("scTypeEval object has correct slots", {
  sceval <- create_test_scTypeEval()
  
  expect_true(methods::is(sceval@counts, "dgCMatrix"))
  expect_true(is.data.frame(sceval@metadata))
  expect_true(is.list(sceval@data))
  expect_true(is.list(sceval@dissimilarity))
  expect_true(is.list(sceval@consistency))
  expect_true(is.list(sceval@gene.lists))
  expect_true(is.list(sceval@reductions))
})


test_that("DataAssay S4 class is defined correctly", {
  sceval <- create_processed_scTypeEval()
  
  expect_s4_class(sceval@data[["single-cell"]], "DataAssay")
})


test_that("DataAssay object has correct slots", {
  sceval <- create_processed_scTypeEval()
  assay <- sceval@data[["single-cell"]]
  
  expect_true(methods::is(assay@matrix, "dgCMatrix"))
  expect_true(is.character(assay@aggregation))
  expect_true(is.factor(assay@group))
  expect_true(is.factor(assay@sample))
  expect_true(is.list(assay@ident))
})


test_that("DissimilarityAssay S4 class is defined correctly", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  expect_s4_class(sceval@dissimilarity[["Pseudobulk:Euclidean"]], "DissimilarityAssay")
})


test_that("DissimilarityAssay object has correct slots", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  assay <- sceval@dissimilarity[["Pseudobulk:Euclidean"]]
  
  expect_s3_class(assay@dissimilarity, "dist")
  expect_true(is.character(assay@method))
  expect_true(is.character(assay@aggregation))
  expect_true(is.factor(assay@sample))
  expect_true(is.list(assay@ident))
})


test_that("DimRed S4 class is defined correctly", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  expect_s4_class(sceval@reductions[["single-cell"]], "DimRed")
})


test_that("DimRed object has correct slots", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  dimred <- sceval@reductions[["single-cell"]]
  
  expect_true(is.matrix(dimred@embeddings))
  expect_true(is.matrix(dimred@feature.loadings))
  expect_true(is.character(dimred@aggregation))
  expect_true(is.factor(dimred@group))
  expect_true(is.factor(dimred@sample))
  expect_true(is.list(dimred@ident))
  expect_true(is.character(dimred@key))
})


test_that("scTypeEval object initialization sets defaults", {
  test_data <- generate_small_test_data()
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  
  expect_equal(length(sceval@data), 0)
  expect_equal(length(sceval@dissimilarity), 0)
  expect_equal(length(sceval@consistency), 0)
  expect_equal(length(sceval@reductions), 0)
  expect_equal(length(sceval@gene.lists), 0)
})


test_that("scTypeEval object stores version information", {
  sceval <- create_test_scTypeEval()
  
  expect_true(is.character(sceval@version))
})


test_that("Multiple aggregation types can coexist", {
  sceval <- create_processed_scTypeEval()
  
  expect_true(all(c("single-cell", "pseudobulk") %in% names(sceval@data)))
})


test_that("Multiple dissimilarity methods can coexist", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Cosine", 
                               reduction = FALSE, verbose = FALSE)
  
  expect_true(all(c("Pseudobulk:Euclidean", "Pseudobulk:Cosine") 
                  %in% names(sceval@dissimilarity)))
})


test_that("Multiple gene lists can coexist", {
  sceval <- create_processed_scTypeEval() # already producing HVG
  sceval <- Run.GeneMarkers(sceval, ngenes.celltype = 50, verbose = FALSE)
  
  expect_true(all(c("HVG", "scran.findMarkers") %in% names(sceval@gene.lists)))
})


test_that("Object structure is preserved after multiple operations", {
  sceval <- create_test_scTypeEval()
  
  # Run multiple operations
  sceval <- set.activeIdent(sceval, ident = "celltype")
  sceval <- Run.ProcessingData(sceval, ident = NULL, sample = "sample", verbose = FALSE)
  sceval <- Run.HVG(sceval, ngenes = 100, verbose = FALSE)
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  # Object should still be valid
  expect_s4_class(sceval, "scTypeEval")
  expect_true(length(sceval@data) > 0)
  expect_true(length(sceval@gene.lists) > 0)
  expect_true(length(sceval@reductions) > 0)
  expect_true(length(sceval@dissimilarity) > 0)
})


test_that("S4 classes maintain data integrity", {
  sceval <- create_processed_scTypeEval()
  
  # Check that dimensions match across related objects
  sc_data <- sceval@data[["single-cell"]]
  expect_equal(length(sc_data@group), ncol(sc_data@matrix))
  expect_equal(length(sc_data@sample), ncol(sc_data@matrix))
  # The ident should have same length as the number of columns (cells)
  expect_equal(length(sc_data@ident[[1]]), ncol(sc_data@matrix))
  # Group and sample should correspond to columns in matrix
  expect_equal(length(sc_data@group), length(sc_data@sample))
})
