test_that("Complete workflow from matrix to consistency metrics", {
  # Create test data with multiple samples
  test_data <- generate_test_data(n_samples = 6)
  
  # Create scTypeEval object
  sceval <- create.scTypeEval(
    matrix = test_data$counts,
    metadata = test_data$metadata
  )
  
  # Process data
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    min.samples = 3,
    min.cells = 5,
    verbose = FALSE
  )
  
  # Identify HVG
  sceval <- Run.HVG(
    sceval,
    var.method = "basic",
    ngenes = 150,
    sample = TRUE,
    verbose = FALSE
  )
  
  # Identify markers
  sceval <- Run.GeneMarkers(
    sceval,
    method = "scran.findMarkers",
    ngenes.celltype = 30,
    verbose = FALSE
  )
  
  # Run PCA
  sceval <- Run.PCA(
    sceval,
    ndim = 5,
    verbose = FALSE
  )
  
  # Compute dissimilarity
  sceval <- Run.Dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = TRUE,
    verbose = FALSE
  )
  
  # Get consistency metrics
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "silhouette",
    verbose = FALSE
  )
  
  # Verify complete workflow
  expect_s4_class(sceval, "scTypeEval")
  expect_true("single-cell" %in% names(sceval@data))
  expect_true("pseudobulk" %in% names(sceval@data))
  expect_true("HVG" %in% names(sceval@gene.lists))
  expect_true("scran.findMarkers" %in% names(sceval@gene.lists))
  expect_true("single-cell" %in% names(sceval@reductions))
  expect_true("pseudobulk" %in% names(sceval@reductions))
  expect_true("Pseudobulk:Euclidean" %in% names(sceval@dissimilarity))
  expect_s3_class(consistency, "data.frame")
  expect_true(nrow(consistency) > 0)
})


test_that("Workflow with different dissimilarity methods", {
  sceval <- create_processed_scTypeEval()
  
  # Test multiple dissimilarity methods
  methods <- c("Pseudobulk:Euclidean", "Pseudobulk:Cosine", "Pseudobulk:Pearson")
  
  for (method in methods) {
    sceval <- Run.Dissimilarity(
      sceval,
      method = method,
      reduction = FALSE,
      verbose = FALSE
    )
  }
  
  # Compute consistency for all methods
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "all",
    Consistency.metric = c("silhouette", "Average.similarity"),
    verbose = FALSE
  )
  
  expect_true(all(methods %in% consistency$dissimilarity_method))
  expect_true(all(c("silhouette", "Average.similarity") %in% consistency$consistency.metric))
})


test_that("Workflow with custom gene lists", {
  sceval <- create_processed_scTypeEval()
  
  # Add custom gene list
  all_genes <- rownames(sceval@data[["single-cell"]]@matrix)
  custom_genes <- all_genes[1:100]
  sceval <- add.GeneList(sceval, gene.list = list("custom" = custom_genes))
  
  # Run PCA with custom genes
  sceval <- Run.PCA(
    sceval,
    gene.list = "custom",
    ndim = 5,
    verbose = FALSE
  )
  
  # Compute dissimilarity
  sceval <- Run.Dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = TRUE,
    verbose = FALSE
  )
  
  expect_equal(sceval@reductions[["pseudobulk"]]@gene.list, custom_genes)
})


test_that("Workflow with black list filtering", {
  sceval <- create_processed_scTypeEval(HVG = FALSE)
  
  # Set black list
  all_genes <- rownames(sceval@data[["single-cell"]]@matrix)
  black_genes <- all_genes[1:20]
  sceval@black.list <- black_genes
  
  # Run HVG with black list
  sceval <- Run.HVG(
    sceval,
    var.method = "basic",
    ngenes = 100,
    verbose = FALSE
  )
  
  # Run PCA with black list
  sceval <- Run.PCA(
    sceval,
    ndim = 5,
    verbose = FALSE
  )
  
  # Verify black list was applied
  expect_true(!any(sceval@gene.lists$HVG %in% black_genes))
})


test_that("Workflow stores consistency in object", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  sceval <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "silhouette",
    return.scTypeEval = TRUE,
    verbose = FALSE
  )
  
  expect_s4_class(sceval, "scTypeEval")
  expect_true(length(sceval@consistency) > 0)
  expect_true("celltype" %in% names(sceval@consistency))
})


test_that("Workflow with Seurat input", {
  skip_if_not_installed("Seurat")
  
  seurat_obj <- create_test_seurat()
  skip_if(is.null(seurat_obj))
  
  # Create from Seurat
  sceval <- create.scTypeEval(seurat_obj)
  
  # Run complete workflow
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    min.samples = 3,
    min.cells = 5,
    verbose = FALSE
  )
  
  sceval <- Run.HVG(sceval, var.method = "basic", ngenes = 100, verbose = FALSE)
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = "silhouette",
    verbose = FALSE
  )
  
  expect_s3_class(consistency, "data.frame")
})


test_that("Workflow with SingleCellExperiment input", {
  skip_if_not_installed("SingleCellExperiment")
  
  sce_obj <- create_test_sce()
  skip_if(is.null(sce_obj))
  
  # Create from SCE
  sceval <- create.scTypeEval(sce_obj)
  
  # Run workflow
  sceval <- Run.ProcessingData(
    sceval,
    ident = "celltype",
    sample = "sample",
    min.samples = 3,
    min.cells = 5,
    verbose = FALSE
  )
  
  sceval <- Run.HVG(sceval, var.method = "basic", ngenes = 100, verbose = FALSE)
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  expect_true("Pseudobulk:Euclidean" %in% names(sceval@dissimilarity))
})


test_that("Workflow with WasserStein dissimilarity", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  
  sceval <- Run.Dissimilarity(
    sceval,
    method = "WasserStein",
    reduction = TRUE,
    verbose = FALSE
  )
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "WasserStein",
    Consistency.metric = c("silhouette", "2label.silhouette"),
    verbose = FALSE
  )
  
  expect_true(all(c("silhouette", "2label.silhouette") %in% consistency$consistency.metric))
})


test_that("Workflow with RecipClassif dissimilarity", {
  skip_if_not_installed("SingleR")
  
  sceval <- create_processed_scTypeEval()
  
  sceval <- Run.Dissimilarity(
    sceval,
    method = "RecipClassif:Match",
    ReciprocalClassifier = "SingleR",
    reduction = FALSE,
    verbose = FALSE
  )
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "RecipClassif:Match",
    Consistency.metric = "silhouette",
    verbose = FALSE
  )
  
  expect_s3_class(consistency, "data.frame")
})


test_that("Workflow preserves metadata throughout", {
  test_data <- generate_test_data(n_samples = 6)
  
  sceval <- create.scTypeEval(test_data$counts, test_data$metadata)
  sceval <- set.activeIdent(sceval, ident = "celltype")
  sceval <- Run.ProcessingData(sceval, ident = NULL, sample = "sample", verbose = FALSE)
  sceval <- Run.HVG(sceval, var.method = "basic", ngenes = 100, verbose = FALSE)
  sceval <- Run.PCA(sceval, ndim = 5, verbose = FALSE)
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  # Check metadata preservation
  expect_equal(nrow(sceval@metadata), ncol(test_data$counts))
  expect_true(all(c("celltype", "sample", "batch", "condition") %in% colnames(sceval@metadata)))
  expect_equal(sceval@active.ident, "celltype")
})


test_that("Complete workflow with all consistency metrics", {
  sceval <- create_processed_scTypeEval()
  sceval <- Run.Dissimilarity(sceval, method = "Pseudobulk:Euclidean", 
                               reduction = FALSE, verbose = FALSE)
  
  all_metrics <- c("silhouette", "2label.silhouette", "NeighborhoodPurity",
                   "ward.PropMatch", "Orbital.medoid", "Average.similarity")
  
  consistency <- get.Consistency(
    sceval,
    dissimilarity.slot = "Pseudobulk:Euclidean",
    Consistency.metric = all_metrics,
    verbose = FALSE
  )
  
  expect_true(all(all_metrics %in% consistency$consistency.metric))
  
  # Each cell type should have measurements for all metrics
  n_celltypes <- length(unique(consistency$celltype))
  expect_equal(nrow(consistency), n_celltypes * length(all_metrics))
})

test_that("Workflow with add.DimReduction integration", {
  sceval <- create_processed_scTypeEval()
  
  # Create custom embeddings
  n_cells <- ncol(sceval@data[["single-cell"]]@matrix)
  custom_emb <- matrix(rnorm(n_cells * 15), nrow = 15, ncol = n_cells)
  
  ident <- sceval@data[["single-cell"]]@ident[[1]]
  sample <- sceval@data[["single-cell"]]@sample
  
  sceval <- add.DimReduction(
    sceval,
    embeddings = custom_emb,
    aggregation = "pseudobulk",
    ident = ident,
    sample = sample,
    key = "custom_reduction",
    verbose = FALSE
  )
  
  # Use custom reduction for dissimilarity
  sceval <- Run.Dissimilarity(
    sceval,
    method = "Pseudobulk:Euclidean",
    reduction = TRUE,
    verbose = FALSE
  )
  
  expect_true("Pseudobulk:Euclidean" %in% names(sceval@dissimilarity))
})
