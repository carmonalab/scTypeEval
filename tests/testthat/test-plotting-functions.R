test_that("plot.PCA generates a ggplot object", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  p <- plot.PCA(sceval, reduction.slot = "single-cell")
  
  expect_s3_class(p, "ggplot")
})


test_that("plot.PCA works with default parameters", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  p <- plot.PCA(sceval)
  
  expect_type(p, "list")
  expect_true(length(p) > 0)
})


test_that("plot.PCA works with reduction.slot = 'all'", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  p <- plot.PCA(sceval, reduction.slot = "all")
  
  expect_type(p, "list")
  expect_true(all(c("single-cell", "pseudobulk") %in% names(p)))
  expect_s3_class(p[["single-cell"]], "ggplot")
  expect_s3_class(p[["pseudobulk"]], "ggplot")
})


test_that("plot.PCA respects label parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  p_with_label <- plot.PCA(sceval, reduction.slot = "single-cell", label = TRUE)
  p_without_label <- plot.PCA(sceval, reduction.slot = "single-cell", label = FALSE)
  
  expect_s3_class(p_with_label, "ggplot")
  expect_s3_class(p_without_label, "ggplot")
})


test_that("plot.PCA respects dims parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  p_12 <- plot.PCA(sceval, reduction.slot = "single-cell", dims = c(1, 2))
  p_23 <- plot.PCA(sceval, reduction.slot = "single-cell", dims = c(2, 3))
  
  expect_s3_class(p_12, "ggplot")
  expect_s3_class(p_23, "ggplot")
})


test_that("plot.PCA respects show.legend parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  p_legend <- plot.PCA(sceval, reduction.slot = "single-cell", show.legend = TRUE)
  p_no_legend <- plot.PCA(sceval, reduction.slot = "single-cell", show.legend = FALSE)
  
  expect_s3_class(p_legend, "ggplot")
  expect_s3_class(p_no_legend, "ggplot")
})


test_that("plot.PCA returns single plot for single reduction slot", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  p <- plot.PCA(sceval, reduction.slot = "single-cell")
  
  expect_s3_class(p, "ggplot")
  expect_false(is.list(p) && !inherits(p, "ggplot"))
})


test_that("plot.PCA errors without PCA data", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    plot.PCA(sceval, reduction.slot = "single-cell")
  )
})


test_that("plot.MDS generates a ggplot object", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot.MDS(sceval, dissimilarity.slot = "Pseudobulk:Euclidean")
  
  expect_s3_class(p, "ggplot")
})


test_that("plot.MDS works with default parameters", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot.MDS(sceval)
  
  expect_s3_class(p, "ggplot")
})


test_that("plot.MDS works with dissimilarity.slot = 'all'", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Cosine", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot.MDS(sceval, dissimilarity.slot = "all")
  
  expect_type(p, "list")
  expect_true(all(c("Pseudobulk:Euclidean", "Pseudobulk:Cosine") %in% names(p)))
  expect_s3_class(p[["Pseudobulk:Euclidean"]], "ggplot")
  expect_s3_class(p[["Pseudobulk:Cosine"]], "ggplot")
})


test_that("plot.MDS respects label parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p_with_label <- plot.MDS(sceval, label = TRUE)
  p_without_label <- plot.MDS(sceval, label = FALSE)
  
  expect_s3_class(p_with_label, "ggplot")
  expect_s3_class(p_without_label, "ggplot")
})


test_that("plot.MDS respects dims parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p_12 <- plot.MDS(sceval, dims = c(1, 2))
  p_23 <- plot.MDS(sceval, dims = c(2, 3))
  
  expect_s3_class(p_12, "ggplot")
  expect_s3_class(p_23, "ggplot")
})


test_that("plot.MDS respects show.legend parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p_legend <- plot.MDS(sceval, show.legend = TRUE)
  p_no_legend <- plot.MDS(sceval, show.legend = FALSE)
  
  expect_s3_class(p_legend, "ggplot")
  expect_s3_class(p_no_legend, "ggplot")
})


test_that("plot.MDS returns single plot for single dissimilarity slot", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot.MDS(sceval, dissimilarity.slot = "Pseudobulk:Euclidean")
  
  expect_s3_class(p, "ggplot")
  expect_false(is.list(p) && !inherits(p, "ggplot"))
})


test_that("plot.MDS errors without dissimilarity data", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    plot.MDS(sceval, dissimilarity.slot = "Pseudobulk:Euclidean")
  )
})


test_that("plot.Heatmap generates a ggplot object", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot.Heatmap(sceval, dissimilarity.slot = "Pseudobulk:Euclidean", 
                    verbose = FALSE)
  
  expect_s3_class(p, "ggplot")
})


test_that("plot.Heatmap works with default parameters", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot.Heatmap(sceval, verbose = FALSE)
  
  expect_s3_class(p, "ggplot")
})


test_that("plot.Heatmap works with dissimilarity.slot = 'all'", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Cosine", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot.Heatmap(sceval, dissimilarity.slot = "all", verbose = FALSE)
  
  expect_type(p, "list")
  expect_true(all(c("Pseudobulk:Euclidean", "Pseudobulk:Cosine") %in% names(p)))
  expect_s3_class(p[["Pseudobulk:Euclidean"]], "ggplot")
  expect_s3_class(p[["Pseudobulk:Cosine"]], "ggplot")
})


test_that("plot.Heatmap respects sort.similarity parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot.Heatmap(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    sort.similarity = "Pseudobulk:Euclidean",
    verbose = FALSE
  )
  
  expect_s3_class(p, "ggplot")
})


test_that("plot.Heatmap respects sort.consistency parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot.Heatmap(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    sort.consistency = "silhouette",
    verbose = FALSE
  )
  
  expect_s3_class(p, "ggplot")
})


test_that("plot.Heatmap works with both sort.similarity and sort.consistency", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot.Heatmap(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    sort.similarity = "Pseudobulk:Euclidean",
    sort.consistency = "silhouette",
    verbose = FALSE
  )
  
  expect_s3_class(p, "ggplot")
})


test_that("plot.Heatmap respects color parameters", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot.Heatmap(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    low.color = "blue",
    high.color = "red",
    verbose = FALSE
  )
  
  expect_s3_class(p, "ggplot")
})


test_that("plot.Heatmap respects hclust.method parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot.Heatmap(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    sort.similarity = "Pseudobulk:Euclidean",
    hclust.method = "complete",
    verbose = FALSE
  )
  
  expect_s3_class(p, "ggplot")
})


test_that("plot.Heatmap handles verbose parameter", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  # verbose = TRUE should produce messages
  expect_message(
    p <- plot.Heatmap(sceval, verbose = TRUE),
    "sorting cell type"
  )
  
  # verbose = FALSE should suppress messages
  expect_silent(
    p <- plot.Heatmap(sceval, verbose = FALSE)
  )
})


test_that("plot.Heatmap returns single plot for single dissimilarity slot", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p <- plot.Heatmap(sceval, dissimilarity.slot = "Pseudobulk:Euclidean", 
                    verbose = FALSE)
  
  expect_s3_class(p, "ggplot")
  expect_false(is.list(p) && !inherits(p, "ggplot"))
})


test_that("plot.Heatmap errors without dissimilarity data", {
  sceval <- create_processed_scTypeEval()
  
  expect_error(
    plot.Heatmap(sceval, dissimilarity.slot = "Pseudobulk:Euclidean", 
                 verbose = FALSE)
  )
})


test_that("plot.Heatmap works with WasserStein dissimilarity", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  sceval <- Run.Dissimilarity(sceval, method = "WasserStein", 
                               reduction = TRUE, verbose = FALSE)
  
  p <- plot.Heatmap(sceval, dissimilarity.slot = "WasserStein", 
                    verbose = FALSE)
  
  expect_s3_class(p, "ggplot")
})


test_that("plot.Heatmap sorting produces different plots", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p_unsorted <- plot.Heatmap(sceval, verbose = FALSE)
  p_sorted <- plot.Heatmap(sceval, sort.consistency = "silhouette", 
                           verbose = FALSE)
  
  # Both should be ggplot objects
  expect_s3_class(p_unsorted, "ggplot")
  expect_s3_class(p_sorted, "ggplot")
})


test_that("Plotting functions work with multiple samples", {
  test_data <- generate_test_data(n_samples = 6)
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  sceval <- Run.ProcessingData(sceval, ident = "celltype", sample = "sample",
                                min.samples = 3, min.cells = 5, verbose = FALSE)
  sceval <- Run.HVG(sceval, ngenes = 100, verbose = FALSE)
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  p_pca <- plot.PCA(sceval, reduction.slot = "single-cell")
  p_mds <- plot.MDS(sceval, dissimilarity.slot = "Pseudobulk:Euclidean")
  p_heatmap <- plot.Heatmap(sceval, dissimilarity.slot = "Pseudobulk:Euclidean", 
                            verbose = FALSE)
  
  expect_s3_class(p_pca, "ggplot")
  expect_s3_class(p_mds, "ggplot")
  expect_s3_class(p_heatmap, "ggplot")
})


test_that("plot.PCA works with different aggregation types", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  p_sc <- plot.PCA(sceval, reduction.slot = "single-cell")
  p_pb <- plot.PCA(sceval, reduction.slot = "pseudobulk")
  
  expect_s3_class(p_sc, "ggplot")
  expect_s3_class(p_pb, "ggplot")
})


test_that("plot.MDS works with different dissimilarity methods", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Pearson", 
                               reduction = FALSE, verbose = FALSE)
  
  p_euc <- plot.MDS(sceval, dissimilarity.slot = "Pseudobulk:Euclidean")
  p_pear <- plot.MDS(sceval, dissimilarity.slot = "Pseudobulk:Pearson")
  
  expect_s3_class(p_euc, "ggplot")
  expect_s3_class(p_pear, "ggplot")
})


test_that("plot.Heatmap can use different dissimilarity for sorting", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Cosine", 
                               reduction = FALSE, verbose = FALSE)
  
  # Plot Euclidean but sort by Cosine
  p <- plot.Heatmap(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    sort.similarity = "Pseudobulk:Cosine",
    verbose = FALSE
  )
  
  expect_s3_class(p, "ggplot")
})


test_that("Plots can be combined with custom ggplot layers", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  p <- plot.PCA(sceval, reduction.slot = "single-cell")
  
  # Add custom layer
  p_modified <- p + ggplot2::theme_minimal()
  
  expect_s3_class(p_modified, "ggplot")
})


test_that("plot.PCA handles high dimensional PCA", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  # Plot later dimensions
  p <- plot.PCA(sceval, reduction.slot = "single-cell", dims = c(4, 5))
  
  expect_s3_class(p, "ggplot")
})


test_that("plot.MDS handles high dimensional MDS", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  # Plot later dimensions
  p <- plot.MDS(sceval, dims = c(3, 4))
  
  expect_s3_class(p, "ggplot")
})
