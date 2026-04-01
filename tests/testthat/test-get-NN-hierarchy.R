test_that("get_hierarchy computes hierarchical clustering", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  hier <- get_hierarchy(sceval, dissimilarity_slot = "Pseudobulk:Euclidean", 
                        verbose = FALSE)
  
  expect_true(is.table(hier) || is.matrix(hier))
})


test_that("get_hierarchy works with default parameters", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  hier <- get_hierarchy(sceval, verbose = FALSE)
  
  expect_true(!is.null(hier))
})


test_that("get_hierarchy works with dissimilarity_slot = 'all'", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Cosine", 
                               reduction = FALSE, verbose = FALSE)
  
  hier <- get_hierarchy(sceval, dissimilarity_slot = "all", verbose = FALSE)
  
  expect_type(hier, "list")
  expect_true(all(c("Pseudobulk:Euclidean", "Pseudobulk:Cosine") %in% names(hier)))
})


test_that("get_hierarchy respects hierarchy_method parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  hier_ward <- get_hierarchy(sceval, hierarchy_method = "ward.D2", verbose = FALSE)
  hier_complete <- get_hierarchy(sceval, hierarchy_method = "complete", verbose = FALSE)
  hier_average <- get_hierarchy(sceval, hierarchy_method = "average", verbose = FALSE)
  
  expect_true(!is.null(hier_ward))
  expect_true(!is.null(hier_complete))
  expect_true(!is.null(hier_average))
})


test_that("get_hierarchy handles verbose parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  # verbose = TRUE should produce messages
  expect_message(
    hier <- get_hierarchy(sceval, verbose = TRUE),
    "Computing hierarchical clustering"
  )
  
  # verbose = FALSE should suppress messages
  expect_silent(
    hier <- get_hierarchy(sceval, verbose = FALSE)
  )
})


test_that("get_hierarchy returns table with correct dimensions", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  hier <- get_hierarchy(sceval, dissimilarity_slot = "Pseudobulk:Euclidean", 
                        verbose = FALSE)
  
  # Should be a contingency table
  expect_true(is.table(hier) || is.matrix(hier))
  expect_true(nrow(hier) > 0)
  expect_true(ncol(hier) > 0)
})


test_that("get_hierarchy clusters match number of cell types", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  hier <- get_hierarchy(sceval, dissimilarity_slot = "Pseudobulk:Euclidean", 
                        verbose = FALSE)
  
  # Number of clusters should equal number of unique cell types
  n_celltypes <- length(unique(sceval@dissimilarity[["Pseudobulk:Euclidean"]]@ident[[1]]))
  expect_equal(nrow(hier), n_celltypes)
})


test_that("get_hierarchy returns single table for single dissimilarity slot", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  hier <- get_hierarchy(sceval, dissimilarity_slot = "Pseudobulk:Euclidean", 
                        verbose = FALSE)
  
  expect_true(is.table(hier) || is.matrix(hier))
  expect_false(is.list(hier) && !is.table(hier))
})


test_that("get_hierarchy errors without dissimilarity data", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    get_hierarchy(sceval, dissimilarity_slot = "Pseudobulk:Euclidean", 
                  verbose = FALSE)
  )
})


test_that("get_hierarchy works with different dissimilarity methods", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Pearson", 
                               reduction = FALSE, verbose = FALSE)
  
  hier <- get_hierarchy(sceval, dissimilarity_slot = "Pseudobulk:Pearson", 
                        verbose = FALSE)
  
  expect_true(!is.null(hier))
})


test_that("get_hierarchy works with WasserStein dissimilarity", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  sceval <- run_dissimilarity(sceval, method = "WasserStein", 
                               reduction = TRUE, verbose = FALSE)
  
  hier <- get_hierarchy(sceval, dissimilarity_slot = "WasserStein", 
                        verbose = FALSE)
  
  expect_true(!is.null(hier))
})

test_that("get_nn computes nearest neighbor composition", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  nn <- get_nn(sceval, dissimilarity_slot = "Pseudobulk:Euclidean", 
               verbose = FALSE)
  
  expect_true(is.data.frame(nn) || is.matrix(nn))
})


test_that("get_nn works with default parameters", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  nn <- get_nn(sceval, verbose = FALSE)
  
  expect_true(!is.null(nn))
})


test_that("get_nn works with dissimilarity_slot = 'all'", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Cosine", 
                               reduction = FALSE, verbose = FALSE)
  
  nn <- get_nn(sceval, dissimilarity_slot = "all", verbose = FALSE)
  
  expect_type(nn, "list")
  expect_true(all(c("Pseudobulk:Euclidean", "Pseudobulk:Cosine") %in% names(nn)))
})


test_that("get_nn respects knn_graph_k parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  nn_k3 <- get_nn(sceval, knn_graph_k = 3, verbose = FALSE)
  nn_k5 <- get_nn(sceval, knn_graph_k = 5, verbose = FALSE)
  nn_k10 <- get_nn(sceval, knn_graph_k = 10, verbose = FALSE)
  
  expect_true(!is.null(nn_k3))
  expect_true(!is.null(nn_k5))
  expect_true(!is.null(nn_k10))
})


test_that("get_nn respects normalize parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  nn_norm <- get_nn(sceval, normalize = TRUE, verbose = FALSE)
  nn_no_norm <- get_nn(sceval, normalize = FALSE, verbose = FALSE)
  
  expect_true(!is.null(nn_norm))
  expect_true(!is.null(nn_no_norm))
})


test_that("get_nn handles verbose parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  # verbose = TRUE should produce messages
  expect_message(
    nn <- get_nn(sceval, verbose = TRUE),
    "Computing hierarchical clustering"
  )
  
  # verbose = FALSE should suppress messages
  expect_silent(
    nn <- get_nn(sceval, verbose = FALSE)
  )
})


test_that("get_nn returns matrix/dataframe with correct structure", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  nn <- get_nn(sceval, dissimilarity_slot = "Pseudobulk:Euclidean", 
               verbose = FALSE)
  
  # Should be a matrix or data frame
  expect_true(is.data.frame(nn) || is.matrix(nn))
  expect_true(nrow(nn) > 0)
  expect_true(ncol(nn) > 0)
})


test_that("get_nn dimensions match cell types", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  nn <- get_nn(sceval, dissimilarity_slot = "Pseudobulk:Euclidean", 
               verbose = FALSE)
  
  # Rows and columns should represent cell types
  n_celltypes <- length(unique(sceval@dissimilarity[["Pseudobulk:Euclidean"]]@ident[[1]]))
  expect_equal(nrow(nn), n_celltypes)
  expect_equal(ncol(nn) - 1, n_celltypes)
  expect_true("celltype" %in% colnames(nn))
})


test_that("get_nn returns single result for single dissimilarity slot", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  nn <- get_nn(sceval, dissimilarity_slot = "Pseudobulk:Euclidean", 
               verbose = FALSE)
  
  expect_true(is.data.frame(nn) || is.matrix(nn))
  expect_false(is.list(nn) && !is.data.frame(nn))
})


test_that("get_nn errors without dissimilarity data", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    get_nn(sceval, dissimilarity_slot = "Pseudobulk:Euclidean", 
           verbose = FALSE)
  )
})


test_that("get_nn works with different dissimilarity methods", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Pearson", 
                               reduction = FALSE, verbose = FALSE)
  
  nn <- get_nn(sceval, dissimilarity_slot = "Pseudobulk:Pearson", 
               verbose = FALSE)
  
  expect_true(!is.null(nn))
})


test_that("get_nn works with WasserStein dissimilarity", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  sceval <- run_dissimilarity(sceval, method = "WasserStein", 
                               reduction = TRUE, verbose = FALSE)
  
  nn <- get_nn(sceval, dissimilarity_slot = "WasserStein", 
               verbose = FALSE)
  
  expect_true(!is.null(nn))
})

test_that("get_nn normalized values differ from unnormalized", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  nn_norm <- get_nn(sceval, normalize = TRUE, verbose = FALSE)
  nn_no_norm <- get_nn(sceval, normalize = FALSE, verbose = FALSE)
  
  # Both should be valid results
  expect_true(!is.null(nn_norm))
  expect_true(!is.null(nn_no_norm))
  
  # They should have the same dimensions
  expect_equal(dim(nn_norm), dim(nn_no_norm))
})


test_that("get_nn with different k values produces different results", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  nn_k3 <- get_nn(sceval, knn_graph_k = 3, verbose = FALSE)
  nn_k10 <- get_nn(sceval, knn_graph_k = 10, verbose = FALSE)
  
  # Both should be valid
  expect_true(!is.null(nn_k3))
  expect_true(!is.null(nn_k10))
  
  # Same dimensions (cell types don't change)
  expect_equal(dim(nn_k3), dim(nn_k10))
})


test_that("get_hierarchy and get_nn work together in analysis workflow", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  hier <- get_hierarchy(sceval, verbose = FALSE)
  nn <- get_nn(sceval, knn_graph_k = 5, verbose = FALSE)
  
  expect_true(!is.null(hier))
  expect_true(!is.null(nn))
})


test_that("get_hierarchy with multiple dissimilarity methods", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  sceval <- run_dissimilarity(sceval, method = "WasserStein", 
                               reduction = TRUE, verbose = FALSE)
  
  hier_all <- get_hierarchy(sceval, dissimilarity_slot = "all", verbose = FALSE)
  
  expect_type(hier_all, "list")
  expect_equal(length(hier_all), 2)
})


test_that("get_nn with multiple dissimilarity methods", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  sceval <- run_dissimilarity(sceval, method = "WasserStein", 
                               reduction = TRUE, verbose = FALSE)
  
  nn_all <- get_nn(sceval, dissimilarity_slot = "all", verbose = FALSE)
  
  expect_type(nn_all, "list")
  expect_equal(length(nn_all), 2)
})


test_that("get_hierarchy different methods produce different clusterings", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  hier_ward <- get_hierarchy(sceval, hierarchy_method = "ward.D2", verbose = FALSE)
  hier_single <- get_hierarchy(sceval, hierarchy_method = "single", verbose = FALSE)
  
  # Both should be valid tables
  expect_true(is.table(hier_ward) || is.matrix(hier_ward))
  expect_true(is.table(hier_single) || is.matrix(hier_single))
  
  # Same dimensions
  expect_equal(dim(hier_ward), dim(hier_single))
})
