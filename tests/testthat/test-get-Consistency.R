test_that("get.Consistency computes silhouette metric", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "silhouette",
    verbose = FALSE
  )
  
  expect_s3_class(consistency, "data.frame")
  expect_true("celltype" %in% colnames(consistency))
  expect_true("measure" %in% colnames(consistency))
  expect_true("consistency.metric" %in% colnames(consistency))
  expect_true("silhouette" %in% consistency$consistency.metric)
})


test_that("get.Consistency computes 2label.silhouette metric", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "2label.silhouette",
    verbose = FALSE
  )
  
  expect_true("2label.silhouette" %in% consistency$consistency.metric)
})


test_that("get.Consistency computes NeighborhoodPurity metric", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "NeighborhoodPurity",
    KNNGraph_k = 5,
    verbose = FALSE
  )
  
  expect_true("NeighborhoodPurity" %in% consistency$consistency.metric)
})


test_that("get.Consistency computes ward.PropMatch metric", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "ward.PropMatch",
    hclust.method = "ward.D2",
    verbose = FALSE
  )
  
  expect_true("ward.PropMatch" %in% consistency$consistency.metric)
})


test_that("get.Consistency computes Orbital.medoid metric", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "Orbital.medoid",
    verbose = FALSE
  )
  
  expect_true("Orbital.medoid" %in% consistency$consistency.metric)
})


test_that("get.Consistency computes Average.similarity metric", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "Average.similarity",
    verbose = FALSE
  )
  
  expect_true("Average.similarity" %in% consistency$consistency.metric)
})


test_that("get.Consistency computes multiple metrics", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = c("silhouette", "2label.silhouette", "Average.similarity"),
    verbose = FALSE
  )
  
  expect_true(all(c("silhouette", "2label.silhouette", "Average.similarity") 
                  %in% consistency$consistency.metric))
})


test_that("get.Consistency works with dissimilarity.slot = 'all'", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Cosine", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "all",
    Consistency.metric = "silhouette",
    verbose = FALSE
  )
  
  expect_true(all(c("Pseudobulk:Euclidean", "Pseudobulk:Cosine") 
                  %in% consistency$dissimilarity_method))
})


test_that("get.Consistency respects KNNGraph_k parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "NeighborhoodPurity",
    KNNGraph_k = 3,
    verbose = FALSE
  )
  
  expect_s3_class(consistency, "data.frame")
})


test_that("get.Consistency respects hclust.method parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "ward.PropMatch",
    hclust.method = "ward.D",
    verbose = FALSE
  )
  
  expect_s3_class(consistency, "data.frame")
})


test_that("get.Consistency respects normalize parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "silhouette",
    normalize = TRUE,
    verbose = FALSE
  )
  
  expect_s3_class(consistency, "data.frame")
})


test_that("get.Consistency handles return.scTypeEval = TRUE", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  sceval_out <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "silhouette",
    return.scTypeEval = TRUE,
    verbose = FALSE
  )
  
  expect_s4_class(sceval_out, "scTypeEval")
  expect_true(length(sceval_out@consistency) > 0)
})


test_that("get.Consistency handles return.scTypeEval = FALSE", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "silhouette",
    return.scTypeEval = FALSE,
    verbose = FALSE
  )
  
  expect_s3_class(consistency, "data.frame")
})


test_that("get.Consistency includes all cell types in results", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "silhouette",
    verbose = FALSE
  )
  
  # Should have results for each cell type
  expect_true(length(unique(consistency$celltype)) > 0)
})


test_that("get.Consistency handles verbose parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  # verbose = TRUE should produce messages
  expect_message(
    consistency <- get.Consistency(
      sceval,
      dissimilarity.slot = "Pseudobulk:Euclidean",
      Consistency.metric = "silhouette",
      verbose = TRUE
    ),
    "Computing internal validation metrics"
  )
  
  # verbose = FALSE should suppress messages
  expect_silent(
    consistency <- get.Consistency(
      sceval,
      dissimilarity.slot = "Pseudobulk:Euclidean",
      Consistency.metric = "silhouette",
      verbose = FALSE
    )
  )
})


test_that("get.Consistency works with multiple samples", {
  test_data <- generate_test_data(n_samples = 6)
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  sceval <- Run.ProcessingData(sceval, ident = "celltype", sample = "sample",
                                min.samples = 3, min.cells = 5, verbose = FALSE)
  sceval <- Run.HVG(sceval, ngenes = 100, verbose = FALSE)
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "silhouette",
    verbose = FALSE
  )
  
  expect_s3_class(consistency, "data.frame")
  expect_true(nrow(consistency) > 0)
})


test_that("get.Consistency returns numeric measures", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "silhouette",
    verbose = FALSE
  )
  
  expect_true(is.numeric(consistency$measure))
})


test_that("get.Consistency errors without dissimilarity data", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    get.Consistency(
      sceval,
      dissimilarity.slot = "Pseudobulk:Euclidean",
      Consistency.metric = "silhouette",
      verbose = FALSE
    )
  )
})
