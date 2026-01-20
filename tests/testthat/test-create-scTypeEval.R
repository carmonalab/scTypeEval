test_that("create.scTypeEval works with count matrix and metadata", {
  test_data <- generate_small_test_data()
  
  sceval <- create.scTypeEval(
    matrix = test_data$counts,
    metadata = test_data$metadata
  )
  
  expect_s4_class(sceval, "scTypeEval")
  expect_equal(ncol(sceval@counts), nrow(sceval@metadata))
  expect_equal(nrow(sceval@counts), nrow(test_data$counts))
  expect_true("celltype" %in% colnames(sceval@metadata))
  expect_true("sample" %in% colnames(sceval@metadata))
})


test_that("create.scTypeEval works with dense matrix", {
  test_data <- generate_small_test_data()
  dense_matrix <- as.matrix(test_data$counts)
  
  sceval <- create.scTypeEval(
    matrix = dense_matrix,
    metadata = test_data$metadata
  )
  
  expect_s4_class(sceval, "scTypeEval")
  expect_s4_class(sceval@counts, "dgCMatrix")
})


test_that("create.scTypeEval works with NULL matrix", {
  metadata <- data.frame(
    celltype = character(0),
    sample = character(0)
  )
  
  expect_message(
    sceval <- create.scTypeEval(matrix = NULL, metadata = metadata),
    "No matrix provided"
  )
  
  expect_s4_class(sceval, "scTypeEval")
  expect_equal(nrow(sceval@counts), 0)
  expect_equal(ncol(sceval@counts), 0)
})


test_that("create.scTypeEval requires metadata for matrix input", {
  test_data <- generate_small_test_data()
  
  expect_error(
    create.scTypeEval(matrix = test_data$counts, metadata = NULL),
    "metadata dataframe must be provided"
  )
})


test_that("create.scTypeEval validates dimensions", {
  test_data <- generate_small_test_data()
  wrong_metadata <- test_data$metadata[1:10, ]
  
  expect_error(
    create.scTypeEval(matrix = test_data$counts, metadata = wrong_metadata),
    "Different number of columns"
  )
})


test_that("create.scTypeEval accepts gene.lists parameter", {
  test_data <- generate_small_test_data()
  gene_list <- list(markers = rownames(test_data$counts)[1:50])
  
  sceval <- create.scTypeEval(
    matrix = test_data$counts,
    metadata = test_data$metadata,
    gene.lists = gene_list
  )
  
  expect_equal(names(sceval@gene.lists), "markers")
  expect_equal(length(sceval@gene.lists$markers), 50)
})


test_that("create.scTypeEval accepts black.list parameter", {
  test_data <- generate_small_test_data()
  black_genes <- rownames(test_data$counts)[1:10]
  
  sceval <- create.scTypeEval(
    matrix = test_data$counts,
    metadata = test_data$metadata,
    black.list = black_genes
  )
  
  expect_equal(length(sceval@black.list), 10)
  expect_true(all(black_genes %in% sceval@black.list))
})


test_that("create.scTypeEval accepts active.ident parameter", {
  test_data <- generate_small_test_data()
  
  sceval <- create.scTypeEval(
    matrix = test_data$counts,
    metadata = test_data$metadata,
    active.ident = "celltype"
  )
  
  expect_equal(sceval@active.ident, "celltype")
})


test_that("create.scTypeEval works with Seurat object", {
  skip_if_not_installed("Seurat")
  
  seurat_obj <- create_test_seurat()
  skip_if(is.null(seurat_obj))
  
  sceval <- create.scTypeEval(seurat_obj)
  
  expect_s4_class(sceval, "scTypeEval")
  expect_s4_class(sceval@counts, "dgCMatrix")
  expect_true(ncol(sceval@counts) > 0)
  expect_true("celltype" %in% colnames(sceval@metadata))
})


test_that("create.scTypeEval works with SingleCellExperiment object", {
  skip_if_not_installed("SingleCellExperiment")
  
  sce_obj <- create_test_sce()
  skip_if(is.null(sce_obj))
  
  sceval <- create.scTypeEval(sce_obj)
  
  expect_s4_class(sceval, "scTypeEval")
  expect_s4_class(sceval@counts, "dgCMatrix")
  expect_true(ncol(sceval@counts) > 0)
  expect_true("celltype" %in% colnames(sceval@metadata))
})


test_that("create.scTypeEval rejects unsupported object types", {
  expect_error(
    create.scTypeEval(matrix = list(a = 1, b = 2)),
    "Input object must be"
  )
})


test_that("create.scTypeEval handles empty gene.lists", {
  test_data <- generate_small_test_data()
  
  sceval <- create.scTypeEval(
    matrix = test_data$counts,
    metadata = test_data$metadata,
    gene.lists = list()
  )
  
  expect_equal(length(sceval@gene.lists), 0)
})


test_that("create.scTypeEval handles multiple samples correctly", {
  test_data <- generate_test_data(n_samples = 6)
  
  sceval <- create.scTypeEval(
    matrix = test_data$counts,
    metadata = test_data$metadata
  )
  
  expect_equal(length(unique(sceval@metadata$sample)), 6)
  expect_true(all(c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6") 
                  %in% unique(sceval@metadata$sample)))
})
