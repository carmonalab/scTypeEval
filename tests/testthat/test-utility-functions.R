test_that("set.activeIdent sets active identity", {
  sceval <- create_test_scTypeEval()
  
  sceval <- set.activeIdent(sceval, ident = "celltype")
  
  expect_equal(sceval@active.ident, "celltype")
})


test_that("set.activeIdent validates ident exists in metadata", {
  sceval <- create_test_scTypeEval()
  
  expect_error(
    set.activeIdent(sceval, ident = "nonexistent_column"),
    "provide a ident"
  )
})


test_that("set.activeIdent requires ident parameter", {
  sceval <- create_test_scTypeEval()
  
  expect_error(
    set.activeIdent(sceval, ident = NULL),
    "Specificy a cell type annotation"
  )
})


test_that("add.GeneList adds a custom gene list", {
  sceval <- create_test_scTypeEval()
  
  all_genes <- rownames(sceval@counts)
  custom_genes <- all_genes[1:50]
  
  sceval <- add.GeneList(
    sceval,
    gene.list = list("custom_markers" = custom_genes)
  )
  
  expect_true("custom_markers" %in% names(sceval@gene.lists))
  expect_equal(length(sceval@gene.lists$custom_markers), 50)
})


test_that("add.GeneList handles multiple gene lists", {
  sceval <- create_test_scTypeEval()
  
  all_genes <- rownames(sceval@counts)
  
  sceval <- add.GeneList(sceval, gene.list = list("list1" = all_genes[1:30]))
  sceval <- add.GeneList(sceval, gene.list = list("list2" = all_genes[31:60]))
  expect_true(all(c("list1", "list2") %in% names(sceval@gene.lists)))
})


test_that("Add.ProcessedData adds single-cell data", {
  sceval <- create_test_scTypeEval()
  test_data <- generate_small_test_data()
  
  # Create test data matching dimensions
  test_counts <- test_data$counts
  test_ident <- test_data$metadata$celltype[1:ncol(test_counts)]
  test_sample <- test_data$metadata$sample[1:ncol(test_counts)]
  
  # Bug fixed: == changed to != in Add.ProcessedData
  # Now correctly-sized data should succeed
  result <- Add.ProcessedData(
    sceval,
    data = test_counts,
    aggregation = "single-cell",
    ident = test_ident,
    ident.name = "celltype",
    sample = test_sample,
    filter = FALSE,
    verbose = FALSE
  )
  
  expect_true("single-cell" %in% names(result@data))
})


test_that("Add.ProcessedData adds pseudobulk data", {
  # Create pseudobulk aggregated data
  test_data <- generate_small_test_data()
  
  # Aggregate by sample and celltype
  pb_data <- test_data$counts
  pb_ident <- test_data$metadata$celltype[1:ncol(pb_data)]
  pb_sample <- test_data$metadata$sample[1:ncol(pb_data)]
  
  sceval <- create_test_scTypeEval()
  # Bug fixed: == changed to != in Add.ProcessedData
  # Now correctly-sized parameters should succeed
  result <- Add.ProcessedData(
    sceval,
    data = pb_data,
    aggregation = "pseudobulk",
    ident = pb_ident,
    sample = pb_sample,
    filter = FALSE,
    verbose = FALSE
  )
  
  expect_true("pseudobulk" %in% names(result@data))
})


test_that("Add.ProcessedData validates ident length", {
  sceval <- create_test_scTypeEval()
  test_data <- generate_small_test_data()
  
  wrong_ident <- rep("CellA", 10)
  
  expect_error(
    Add.ProcessedData(
      sceval,
      data = test_data$counts,
      aggregation = "single-cell",
      ident = wrong_ident,
      sample = test_data$metadata$sample,
      verbose = FALSE
    )
  )
})


test_that("Add.ProcessedData validates sample length", {
  sceval <- create_test_scTypeEval()
  test_data <- generate_small_test_data()
  
  wrong_sample <- rep("Sample1", 10)
  
  expect_error(
    Add.ProcessedData(
      sceval,
      data = test_data$counts,
      aggregation = "single-cell",
      ident = test_data$metadata$celltype,
      sample = wrong_sample,
      verbose = FALSE
    )
  )
})


test_that("Add.ProcessedData handles filter parameter", {
  sceval <- create_test_scTypeEval()
  test_data <- generate_small_test_data()
  
  # Bug fixed: == changed to != in Add.ProcessedData
  # Now correctly-sized parameters with filter should succeed
  result <- Add.ProcessedData(
    sceval,
    data = test_data$counts,
    aggregation = "single-cell",
    ident = test_data$metadata$celltype,
    sample = test_data$metadata$sample,
    filter = TRUE,
    min.samples = 2,
    min.cells = 5,
    verbose = FALSE
  )
  
  expect_true("single-cell" %in% names(result@data))
})


test_that("add.DimReduction adds custom dimensional reduction", {
  sceval <- create_processed_scTypeEval()
  
  # Create mock embeddings
  n_cells <- ncol(sceval@data[["single-cell"]]@matrix)
  embeddings <- matrix(rnorm(n_cells * 10), nrow = 10, ncol = n_cells)
  
  ident <- sceval@data[["single-cell"]]@ident[[1]]
  sample <- sceval@data[["single-cell"]]@sample
  
  sceval <- add.DimReduction(
    sceval,
    embeddings = embeddings,
    aggregation = "single-cell",
    ident = ident,
    sample = sample,
    key = "custom_embedding",
    verbose = FALSE
  )
  
  expect_true("single-cell" %in% names(sceval@reductions))
  expect_equal(sceval@reductions[["single-cell"]]@key, "custom_embedding")
})


test_that("add.DimReduction validates embeddings dimensions", {
  sceval <- create_processed_scTypeEval()
  
  # Create wrong-sized embeddings
  embeddings <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  
  ident <- sceval@data[["single-cell"]]@ident[[1]]
  sample <- sceval@data[["single-cell"]]@sample
  
  expect_error(
    add.DimReduction(
      sceval,
      embeddings = embeddings,
      aggregation = "single-cell",
      ident = ident,
      sample = sample,
      key = "custom",
      verbose = FALSE
    )
  )
})


test_that("add.DimReduction accepts feature.loadings parameter", {
  sceval <- create_processed_scTypeEval()
  
  n_cells <- ncol(sceval@data[["single-cell"]]@matrix)
  n_genes <- nrow(sceval@data[["single-cell"]]@matrix)
  embeddings <- matrix(rnorm(n_cells * 10), nrow = 10, ncol = n_cells)
  loadings <- matrix(rnorm(n_genes * 10), nrow = n_genes, ncol = 10)
  
  ident <- sceval@data[["single-cell"]]@ident[[1]]
  sample <- sceval@data[["single-cell"]]@sample
  
  sceval <- add.DimReduction(
    sceval,
    embeddings = embeddings,
    aggregation = "single-cell",
    ident = ident,
    sample = sample,
    key = "custom_pca",
    feature.loadings = loadings,
    verbose = FALSE
  )
  
  expect_true(nrow(sceval@reductions[["single-cell"]]@feature.loadings) > 0)
})


test_that("add.DimReduction handles gene.list parameter", {
  sceval <- create_processed_scTypeEval() # already producing HVG
  
  n_cells <- ncol(sceval@data[["single-cell"]]@matrix)
  embeddings <- matrix(rnorm(n_cells * 10), nrow = 10, ncol = n_cells)
  
  ident <- sceval@data[["single-cell"]]@ident[[1]]
  sample <- sceval@data[["single-cell"]]@sample
  
  sceval <- add.DimReduction(
    sceval,
    embeddings = embeddings,
    aggregation = "single-cell",
    ident = ident,
    sample = sample,
    key = "custom",
    gene.list = "HVG",
    verbose = FALSE
  )
  
  expect_equal(sceval@reductions[["single-cell"]]@gene.list, "HVG")
})


test_that("add.DimReduction handles black.list parameter", {
  sceval <- create_processed_scTypeEval()
  
  all_genes <- rownames(sceval@data[["single-cell"]]@matrix)
  black_genes <- all_genes[1:10]
  
  n_cells <- ncol(sceval@data[["single-cell"]]@matrix)
  embeddings <- matrix(rnorm(n_cells * 10), nrow = 10, ncol = n_cells)
  
  ident <- sceval@data[["single-cell"]]@ident[[1]]
  sample <- sceval@data[["single-cell"]]@sample
  
  sceval <- add.DimReduction(
    sceval,
    embeddings = embeddings,
    aggregation = "single-cell",
    ident = ident,
    sample = sample,
    key = "custom",
    black.list = black_genes,
    verbose = FALSE
  )
  
  expect_true(length(sceval@reductions[["single-cell"]]@black.list) > 0)
})


test_that("utility functions work together", {
  sceval <- create_test_scTypeEval()
  
  # Set active ident
  sceval <- set.activeIdent(sceval, ident = "celltype")
  
  # Add custom gene list (use correct format: list with named elements)
  all_genes <- rownames(sceval@counts)
  sceval <- add.GeneList(sceval, gene.list = list("custom" = all_genes[1:50]))
  
  # Process data using active ident
  sceval <- Run.ProcessingData(
    sceval,
    ident = NULL,  # Should use active.ident
    sample = "sample",
    verbose = FALSE
  )
  
  expect_true("single-cell" %in% names(sceval@data))
  expect_true("custom" %in% names(sceval@gene.lists))
  expect_equal(sceval@active.ident, "celltype")
})
