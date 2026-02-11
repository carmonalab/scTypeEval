# Test suite for load_singleCell_object function
# 
# This file tests the load_singleCell_object function which handles
# loading single-cell data from various formats (.rds, .h5ad)

# Helper function to create test files ---------------------------------------

create_test_seurat_file <- function(filename) {
  skip_if_not_installed("Seurat")
  
  test_dir <- file.path(tempdir(), "scTypeEval_load_tests")
  dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)
  filepath <- file.path(test_dir, filename)
  
  # Create small Seurat object with sparse matrix
  counts <- Matrix::Matrix(rpois(1000, 5), nrow = 50, ncol = 20, sparse = TRUE)
  rownames(counts) <- paste0("Gene", 1:50)
  colnames(counts) <- paste0("Cell", 1:20)
  
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = data.frame(
      cell_type = rep(c("TypeA", "TypeB"), each = 10),
      sample = rep(paste0("Sample", 1:4), each = 5),
      row.names = colnames(counts)
    )
  )
  
  saveRDS(seurat_obj, filepath)
  return(filepath)
}

create_test_sce_file <- function(filename) {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  test_dir <- file.path(tempdir(), "scTypeEval_load_tests")
  dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)
  filepath <- file.path(test_dir, filename)
  
  # Create small SCE object with sparse matrix
  counts <- Matrix::Matrix(rpois(1000, 5), nrow = 50, ncol = 20, sparse = TRUE)
  rownames(counts) <- paste0("Gene", 1:50)
  colnames(counts) <- paste0("Cell", 1:20)
  
  sce_obj <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = data.frame(
      cell_type = rep(c("TypeA", "TypeB"), each = 10),
      sample = rep(paste0("Sample", 1:4), each = 5),
      row.names = colnames(counts)
    )
  )
  
  saveRDS(sce_obj, filepath)
  return(filepath)
}

# Tests for Seurat objects ---------------------------------------------------

test_that("load_singleCell_object loads Seurat object with split=TRUE", {
  skip_if_not_installed("Seurat")
  
  filepath <- create_test_seurat_file("seurat_test.rds")
  
  result <- load_singleCell_object(filepath, split = TRUE)
  
  expect_type(result, "list")
  expect_named(result, c("counts", "metadata"))
  expect_s4_class(result$counts, "dgCMatrix")
  expect_s3_class(result$metadata, "data.frame")
  expect_equal(ncol(result$counts), nrow(result$metadata))
  expect_equal(ncol(result$counts), 20)
  expect_equal(nrow(result$counts), 50)
})

test_that("load_singleCell_object loads Seurat object with split=FALSE", {
  skip_if_not_installed("Seurat")
  
  filepath <- create_test_seurat_file("seurat_test2.rds")
  
  result <- load_singleCell_object(filepath, split = FALSE)
  
  expect_s4_class(result, "Seurat")
  expect_equal(ncol(result), 20)
  expect_equal(nrow(result), 50)
})

test_that("Seurat object returns correct metadata columns", {
  skip_if_not_installed("Seurat")
  
  filepath <- create_test_seurat_file("seurat_test3.rds")
  
  result <- load_singleCell_object(filepath, split = TRUE)
  
  expect_true("cell_type" %in% colnames(result$metadata))
  expect_true("sample" %in% colnames(result$metadata))
  expect_equal(length(unique(result$metadata$cell_type)), 2)
  expect_equal(length(unique(result$metadata$sample)), 4)
})

test_that("Seurat counts matrix has correct format", {
  skip_if_not_installed("Seurat")
  
  filepath <- create_test_seurat_file("seurat_test4.rds")
  
  result <- load_singleCell_object(filepath, split = TRUE)
  
  expect_true(inherits(result$counts, "dgCMatrix"))
  expect_true(all(result$counts >= 0))
  expect_equal(colnames(result$counts), rownames(result$metadata))
})

# Tests for SingleCellExperiment objects -------------------------------------

test_that("load_singleCell_object loads SCE object with split=TRUE", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  filepath <- create_test_sce_file("sce_test.rds")
  
  result <- load_singleCell_object(filepath, split = TRUE)
  
  expect_type(result, "list")
  expect_named(result, c("counts", "metadata"))
  expect_s4_class(result$counts, "dgCMatrix")
  expect_s3_class(result$metadata, "data.frame")
  expect_equal(ncol(result$counts), nrow(result$metadata))
  expect_equal(ncol(result$counts), 20)
  expect_equal(nrow(result$counts), 50)
})

test_that("load_singleCell_object loads SCE object with split=FALSE", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  filepath <- create_test_sce_file("sce_test2.rds")
  
  result <- load_singleCell_object(filepath, split = FALSE)
  
  expect_s4_class(result, "SingleCellExperiment")
  expect_equal(ncol(result), 20)
  expect_equal(nrow(result), 50)
})

test_that("SCE object returns correct metadata columns", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  filepath <- create_test_sce_file("sce_test3.rds")
  
  result <- load_singleCell_object(filepath, split = TRUE)
  
  expect_true("cell_type" %in% colnames(result$metadata))
  expect_true("sample" %in% colnames(result$metadata))
  expect_equal(length(unique(result$metadata$cell_type)), 2)
  expect_equal(length(unique(result$metadata$sample)), 4)
})

test_that("SCE counts matrix has correct format", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  filepath <- create_test_sce_file("sce_test4.rds")
  
  result <- load_singleCell_object(filepath, split = TRUE)
  
  expect_true(inherits(result$counts, "dgCMatrix"))
  expect_true(all(result$counts >= 0))
  expect_equal(colnames(result$counts), rownames(result$metadata))
})

# Tests for error handling ---------------------------------------------------

test_that("load_singleCell_object errors when file does not exist", {
  nonexistent_file <- file.path(tempdir(), "nonexistent_file.rds")
  
  expect_error(
    load_singleCell_object(nonexistent_file, split = TRUE),
    "File does not exist"
  )
})

test_that("load_singleCell_object errors with unsupported .rds object type", {
  test_dir <- file.path(tempdir(), "scTypeEval_load_tests")
  dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)
  filepath <- file.path(test_dir, "unsupported.rds")
  
  # Save a regular list (unsupported type)
  saveRDS(list(data = 1:10), filepath)
  
  expect_error(
    load_singleCell_object(filepath, split = TRUE),
    "Unsupported .rds object type"
  )
})

test_that("load_singleCell_object errors with unsupported file format", {
  test_dir <- file.path(tempdir(), "scTypeEval_load_tests")
  dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)
  filepath <- file.path(test_dir, "test.txt")
  
  # Create a text file
  writeLines("test", filepath)
  
  expect_error(
    load_singleCell_object(filepath, split = TRUE),
    "Unsupported file format"
  )
})

# Tests for parameter validation ---------------------------------------------

test_that("load_singleCell_object accepts boolean split parameter", {
  skip_if_not_installed("Seurat")
  
  filepath <- create_test_seurat_file("seurat_test_bool.rds")
  
  # Test with TRUE
  result_true <- load_singleCell_object(filepath, split = TRUE)
  expect_type(result_true, "list")
  
  # Test with FALSE
  result_false <- load_singleCell_object(filepath, split = FALSE)
  expect_s4_class(result_false, "Seurat")
})

# Integration tests ----------------------------------------------------------

test_that("loaded Seurat data can create scTypeEval object", {
  skip_if_not_installed("Seurat")
  
  filepath <- create_test_seurat_file("seurat_integration.rds")
  
  result <- load_singleCell_object(filepath, split = TRUE)
  
  # Create scTypeEval object from loaded data
  obj <- create.scTypeEval(
    matrix = result$counts,
    metadata = result$metadata,
    active.ident = "cell_type"
  )
  obj <- set.activeIdent(obj, ident = "sample")
  
  expect_s4_class(obj, "scTypeEval")
  expect_equal(ncol(obj@counts), 20)
})

test_that("loaded SCE data can create scTypeEval object", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  filepath <- create_test_sce_file("sce_integration.rds")
  
  result <- load_singleCell_object(filepath, split = TRUE)
  
  # Create scTypeEval object from loaded data
  obj <- create.scTypeEval(
    matrix = result$counts,
    metadata = result$metadata,
    active.ident = "cell_type"
  )
  obj <- set.activeIdent(obj, ident = "sample")
  
  expect_s4_class(obj, "scTypeEval")
  expect_equal(ncol(obj@counts), 20)
})

test_that("loaded Seurat object directly can create scTypeEval object", {
  skip_if_not_installed("Seurat")
  
  filepath <- create_test_seurat_file("seurat_direct.rds")
  
  seurat_obj <- load_singleCell_object(filepath, split = FALSE)
  
  # Create scTypeEval object from Seurat object directly
  obj <- create.scTypeEval(
    matrix = seurat_obj,
    active.ident = "cell_type"
  )
  obj <- set.activeIdent(obj, ident = "sample")
  
  expect_s4_class(obj, "scTypeEval")
  expect_equal(ncol(obj@counts), 20)
})

test_that("loaded SCE object directly can create scTypeEval object", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  filepath <- create_test_sce_file("sce_direct.rds")
  
  sce_obj <- load_singleCell_object(filepath, split = FALSE)
  
  # Create scTypeEval object from SCE object directly
  obj <- create.scTypeEval(
    matrix = sce_obj,
    active.ident = "cell_type"
  )
  obj <- set.activeIdent(obj, ident = "sample")
  
  expect_s4_class(obj, "scTypeEval")
  expect_equal(ncol(obj@counts), 20)
})

# Tests for data integrity ---------------------------------------------------

test_that("Seurat loaded data preserves gene and cell names", {
  skip_if_not_installed("Seurat")
  
  filepath <- create_test_seurat_file("seurat_names.rds")
  
  result <- load_singleCell_object(filepath, split = TRUE)
  
  expect_true(all(grepl("^Gene", rownames(result$counts))))
  expect_true(all(grepl("^Cell", colnames(result$counts))))
  expect_true(all(grepl("^Cell", rownames(result$metadata))))
})

test_that("SCE loaded data preserves gene and cell names", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  
  filepath <- create_test_sce_file("sce_names.rds")
  
  result <- load_singleCell_object(filepath, split = TRUE)
  
  expect_true(all(grepl("^Gene", rownames(result$counts))))
  expect_true(all(grepl("^Cell", colnames(result$counts))))
  expect_true(all(grepl("^Cell", rownames(result$metadata))))
})

test_that("Loaded data has matching dimensions between counts and metadata", {
  skip_if_not_installed("Seurat")
  
  filepath <- create_test_seurat_file("seurat_dims.rds")
  
  result <- load_singleCell_object(filepath, split = TRUE)
  
  expect_equal(ncol(result$counts), nrow(result$metadata))
  expect_equal(colnames(result$counts), rownames(result$metadata))
})

test_that("Loaded counts matrix contains only non-negative values", {
  skip_if_not_installed("Seurat")
  
  filepath <- create_test_seurat_file("seurat_values.rds")
  
  result <- load_singleCell_object(filepath, split = TRUE)
  
  expect_true(all(result$counts@x >= 0))
  expect_true(min(result$counts) >= 0)
})
