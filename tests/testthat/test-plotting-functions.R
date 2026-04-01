test_that("plot_pca generates a ggplot object", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  p <- plot_pca(sceval, reduction_slot = "single-cell")
  
  expect_s3_class(p, "ggplot")
})


test_that("plot_pca works with default parameters", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  p <- plot_pca(sceval)
  
  expect_type(p, "list")
  expect_true(length(p) > 0)
})


test_that("plot_pca works with reduction_slot = 'all'", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  p <- plot_pca(sceval, reduction_slot = "all")
  
  expect_type(p, "list")
  expect_true(all(c("single-cell", "pseudobulk") %in% names(p)))
  expect_s3_class(p[["single-cell"]], "ggplot")
  expect_s3_class(p[["pseudobulk"]], "ggplot")
})


test_that("plot_pca respects label parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  p_with_label <- plot_pca(sceval, reduction_slot = "single-cell", label = TRUE)
  p_without_label <- plot_pca(sceval, reduction_slot = "single-cell", label = FALSE)
  
  expect_s3_class(p_with_label, "ggplot")
  expect_s3_class(p_without_label, "ggplot")
})


test_that("plot_pca respects dims parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  p_12 <- plot_pca(sceval, reduction_slot = "single-cell", dims = c(1, 2))
  p_23 <- plot_pca(sceval, reduction_slot = "single-cell", dims = c(2, 3))
  
  expect_s3_class(p_12, "ggplot")
  expect_s3_class(p_23, "ggplot")
})


test_that("plot_pca respects show_legend parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  p_legend <- plot_pca(sceval, reduction_slot = "single-cell", show_legend = TRUE)
  p_no_legend <- plot_pca(sceval, reduction_slot = "single-cell", show_legend = FALSE)
  
  expect_s3_class(p_legend, "ggplot")
  expect_s3_class(p_no_legend, "ggplot")
})


test_that("plot_pca returns single plot for single reduction slot", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  p <- plot_pca(sceval, reduction_slot = "single-cell")
  
  expect_s3_class(p, "ggplot")
  expect_false(is.list(p) && !inherits(p, "ggplot"))
})


test_that("plot_pca errors without PCA data", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    plot_pca(sceval, reduction_slot = "single-cell")
  )
})


test_that("plot_mds generates a ggplot object", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot_mds(sceval, dissimilarity_slot = "Pseudobulk:Euclidean")
  
  expect_s3_class(p, "ggplot")
})


test_that("plot_mds works with default parameters", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot_mds(sceval)
  
  expect_s3_class(p, "ggplot")
})


test_that("plot_mds works with dissimilarity_slot = 'all'", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Cosine", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot_mds(sceval, dissimilarity_slot = "all")
  
  expect_type(p, "list")
  expect_true(all(c("Pseudobulk:Euclidean", "Pseudobulk:Cosine") %in% names(p)))
  expect_s3_class(p[["Pseudobulk:Euclidean"]], "ggplot")
  expect_s3_class(p[["Pseudobulk:Cosine"]], "ggplot")
})


test_that("plot_mds respects label parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p_with_label <- plot_mds(sceval, label = TRUE)
  p_without_label <- plot_mds(sceval, label = FALSE)
  
  expect_s3_class(p_with_label, "ggplot")
  expect_s3_class(p_without_label, "ggplot")
})


test_that("plot_mds respects dims parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p_12 <- plot_mds(sceval, dims = c(1, 2))
  p_23 <- plot_mds(sceval, dims = c(2, 3))
  
  expect_s3_class(p_12, "ggplot")
  expect_s3_class(p_23, "ggplot")
})


test_that("plot_mds respects show_legend parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p_legend <- plot_mds(sceval, show_legend = TRUE)
  p_no_legend <- plot_mds(sceval, show_legend = FALSE)
  
  expect_s3_class(p_legend, "ggplot")
  expect_s3_class(p_no_legend, "ggplot")
})


test_that("plot_mds returns single plot for single dissimilarity slot", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot_mds(sceval, dissimilarity_slot = "Pseudobulk:Euclidean")
  
  expect_s3_class(p, "ggplot")
  expect_false(is.list(p) && !inherits(p, "ggplot"))
})


test_that("plot_mds errors without dissimilarity data", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    plot_mds(sceval, dissimilarity_slot = "Pseudobulk:Euclidean")
  )
})


test_that("plot_heatmap generates a ggplot object", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot_heatmap(sceval, dissimilarity_slot = "Pseudobulk:Euclidean", 
                    verbose = FALSE)
  
  expect_s3_class(p, "ggplot")
})


test_that("plot_heatmap works with default parameters", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot_heatmap(sceval, verbose = FALSE)
  
  expect_s3_class(p, "ggplot")
})


test_that("plot_heatmap works with dissimilarity_slot = 'all'", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Cosine", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot_heatmap(sceval, dissimilarity_slot = "all", verbose = FALSE)
  
  expect_type(p, "list")
  expect_true(all(c("Pseudobulk:Euclidean", "Pseudobulk:Cosine") %in% names(p)))
  expect_s3_class(p[["Pseudobulk:Euclidean"]], "ggplot")
  expect_s3_class(p[["Pseudobulk:Cosine"]], "ggplot")
})


test_that("plot_heatmap respects sort_similarity parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot_heatmap(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    sort_similarity = "Pseudobulk:Euclidean",
    verbose = FALSE
  )
  
  expect_s3_class(p, "ggplot")
})


test_that("plot_heatmap respects sort_consistency parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot_heatmap(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    sort_consistency = "silhouette",
    verbose = FALSE
  )
  
  expect_s3_class(p, "ggplot")
})


test_that("plot_heatmap works with both sort_similarity and sort_consistency", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot_heatmap(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    sort_similarity = "Pseudobulk:Euclidean",
    sort_consistency = "silhouette",
    verbose = FALSE
  )
  
  expect_s3_class(p, "ggplot")
})


test_that("plot_heatmap respects color parameters", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot_heatmap(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    low_color = "blue",
    high_color = "red",
    verbose = FALSE
  )
  
  expect_s3_class(p, "ggplot")
})


test_that("plot_heatmap respects hclust_method parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot_heatmap(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    sort_similarity = "Pseudobulk:Euclidean",
    hclust_method = "complete",
    verbose = FALSE
  )
  
  expect_s3_class(p, "ggplot")
})


test_that("plot_heatmap handles verbose parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  # verbose = TRUE should produce messages
  expect_message(
    p <- plot_heatmap(sceval, verbose = TRUE),
    "sorting cell type"
  )
  
  # verbose = FALSE should suppress messages
  expect_silent(
    p <- plot_heatmap(sceval, verbose = FALSE)
  )
})


test_that("plot_heatmap returns single plot for single dissimilarity slot", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot_heatmap(sceval, dissimilarity_slot = "Pseudobulk:Euclidean", 
                    verbose = FALSE)
  
  expect_s3_class(p, "ggplot")
  expect_false(is.list(p) && !inherits(p, "ggplot"))
})


test_that("plot_heatmap errors without dissimilarity data", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    plot_heatmap(sceval, dissimilarity_slot = "Pseudobulk:Euclidean", 
                 verbose = FALSE)
  )
})


test_that("plot_heatmap works with WasserStein dissimilarity", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  sceval <- run_dissimilarity(sceval, method = "WasserStein", 
                               reduction = TRUE, verbose = FALSE)
  
  p <- plot_heatmap(sceval, dissimilarity_slot = "WasserStein", 
                    verbose = FALSE)
  
  expect_s3_class(p, "ggplot")
})


test_that("plot_heatmap sorting produces different plots", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p_unsorted <- plot_heatmap(sceval, verbose = FALSE)
  p_sorted <- plot_heatmap(sceval, sort_consistency = "silhouette", 
                           verbose = FALSE)
  
  # Both should be ggplot objects
  expect_s3_class(p_unsorted, "ggplot")
  expect_s3_class(p_sorted, "ggplot")
})

test_that("plot_pca works with different aggregation types", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  p_sc <- plot_pca(sceval, reduction_slot = "single-cell")
  p_pb <- plot_pca(sceval, reduction_slot = "pseudobulk")
  
  expect_s3_class(p_sc, "ggplot")
  expect_s3_class(p_pb, "ggplot")
})


test_that("plot_mds works with different dissimilarity methods", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Pearson", 
                               reduction = FALSE, verbose = FALSE)
  
  p_euc <- plot_mds(sceval, dissimilarity_slot = "Pseudobulk:Euclidean")
  p_pear <- plot_mds(sceval, dissimilarity_slot = "Pseudobulk:Pearson")
  
  expect_s3_class(p_euc, "ggplot")
  expect_s3_class(p_pear, "ggplot")
})


test_that("plot_heatmap can use different dissimilarity for sorting", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Cosine", 
                               reduction = FALSE, verbose = FALSE)
  
  # Plot Euclidean but sort by Cosine
  p <- plot_heatmap(
    sceval,
    dissimilarity_slot = "Pseudobulk:Euclidean",
    sort_similarity = "Pseudobulk:Cosine",
    verbose = FALSE
  )
  
  expect_s3_class(p, "ggplot")
})


test_that("Plots can be combined with custom ggplot layers", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  p <- plot_pca(sceval, reduction_slot = "single-cell")
  
  # Add custom layer
  p_modified <- p + ggplot2::theme_minimal()
  
  expect_s3_class(p_modified, "ggplot")
})


test_that("plot_pca handles high dimensional PCA", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_pca(sceval, ndim = 5, verbose = FALSE)
  
  # Plot later dimensions
  p <- plot_pca(sceval, reduction_slot = "single-cell", dims = c(4, 5))
  
  expect_s3_class(p, "ggplot")
})


test_that("plot_mds handles high dimensional MDS", {
  sceval <- create_processed_scTypeEval()
  sceval <- run_dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  # Plot later dimensions
  p <- plot_mds(sceval, dims = c(3, 4))
  
  expect_s3_class(p, "ggplot")
})
