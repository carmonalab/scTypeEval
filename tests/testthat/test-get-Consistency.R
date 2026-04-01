test_that("get_consistency computes silhouette metric", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get_consistency(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    consistency_metric = "silhouette",
    verbose = FALSE
  )
  
  expect_s3_class(consistency, "data.frame")
  expect_true("celltype" %in% colnames(consistency))
  expect_true("measure" %in% colnames(consistency))
  expect_true("consistency_metric" %in% colnames(consistency))
  expect_true("silhouette" %in% consistency$consistency_metric)
})


test_that("get_consistency computes 2label_silhouette metric", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get_consistency(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    consistency_metric = "2label_silhouette",
    verbose = FALSE
  )
  
  expect_true("2label_silhouette" %in% consistency$consistency_metric)
})


test_that("get_consistency computes NeighborhoodPurity metric", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get_consistency(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    consistency_metric = "NeighborhoodPurity",
    knn_graph_k = 5,
    verbose = FALSE
  )
  
  expect_true("NeighborhoodPurity" %in% consistency$consistency_metric)
})


test_that("get_consistency computes ward_PropMatch metric", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get_consistency(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    consistency_metric = "ward_PropMatch",
    hclust_method = "ward.D2",
    verbose = FALSE
  )
  
  expect_true("ward_PropMatch" %in% consistency$consistency_metric)
})


test_that("get_consistency computes Orbital_medoid metric", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get_consistency(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    consistency_metric = "Orbital_medoid",
    verbose = FALSE
  )
  
  expect_true("Orbital_medoid" %in% consistency$consistency_metric)
})


test_that("get_consistency computes Average_similarity metric", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get_consistency(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    consistency_metric = "Average_similarity",
    verbose = FALSE
  )
  
  expect_true("Average_similarity" %in% consistency$consistency_metric)
})


test_that("get_consistency computes multiple metrics", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get_consistency(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    consistency_metric = c("silhouette", "2label_silhouette", "Average_similarity"),
    verbose = FALSE
  )
  
  expect_true(all(c("silhouette", "2label_silhouette", "Average_similarity") 
                   %in% consistency$consistency_metric))
})


test_that("get_consistency works with dissimilarity_slot = 'all'", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Cosine", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get_consistency(
    sceval,
    dissimilarity_slot = "all",
    consistency_metric = "silhouette",
    verbose = FALSE
  )
  
  expect_true(all(c("Pseudobulk:Euclidean", "Pseudobulk:Cosine") 
                  %in% consistency$dissimilarity_method))
})


test_that("get_consistency respects knn_graph_k parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get_consistency(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    consistency_metric = "NeighborhoodPurity",
    knn_graph_k = 3,
    verbose = FALSE
  )
  
  expect_s3_class(consistency, "data.frame")
})


test_that("get_consistency respects hclust_method parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get_consistency(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    consistency_metric = "ward_PropMatch",
    hclust_method = "ward.D",
    verbose = FALSE
  )
  
  expect_s3_class(consistency, "data.frame")
})


test_that("get_consistency respects normalize parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get_consistency(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    consistency_metric = "silhouette",
    normalize = TRUE,
    verbose = FALSE
  )
  
  expect_s3_class(consistency, "data.frame")
})


test_that("get_consistency handles return_scTypeEval = TRUE", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  sceval_out <- get_consistency(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    consistency_metric = "silhouette",
    return_scTypeEval = TRUE,
    verbose = FALSE
  )
  
  expect_s4_class(sceval_out, "scTypeEval")
  expect_true(length(sceval_out@consistency) > 0)
})


test_that("get_consistency handles return_scTypeEval = FALSE", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get_consistency(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    consistency_metric = "silhouette",
    return_scTypeEval = FALSE,
    verbose = FALSE
  )
  
  expect_s3_class(consistency, "data.frame")
})


test_that("get_consistency includes all cell types in results", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get_consistency(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    consistency_metric = "silhouette",
    verbose = FALSE
  )
  
  # Should have results for each cell type
  expect_true(length(unique(consistency$celltype)) > 0)
})


test_that("get_consistency handles verbose parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  # verbose = TRUE should produce messages
  expect_message(
    consistency <- get_consistency(
      sceval,
      dissimilarity_slot = "Pseudobulk:Euclidean",
      consistency_metric = "silhouette",
      verbose = TRUE
    ),
    "Computing internal validation metrics"
  )
  
  # verbose = FALSE should suppress messages
  expect_silent(
    consistency <- get_consistency(
      sceval,
      dissimilarity_slot = "Pseudobulk:Euclidean",
      consistency_metric = "silhouette",
      verbose = FALSE
    )
  )
})


test_that("get_consistency returns numeric measures", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get_consistency(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    consistency_metric = "silhouette",
    verbose = FALSE
  )
  
  expect_true(is.numeric(consistency$measure))
})


test_that("get_consistency errors without dissimilarity data", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    get_consistency(
      sceval,
      dissimilarity_slot = "Pseudobulk:Euclidean",
      consistency_metric = "silhouette",
      verbose = FALSE
    )
  )
})
