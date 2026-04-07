# scTypeEval News

## Version 0.99.30 (2026-04-07)

### Improvements

* Updated both vignettes to use fully runnable toy data workflows for matrix, Seurat, and SingleCellExperiment inputs.
* Removed non-justified `eval = FALSE` vignette chunks and aligned examples with the current API.
* Updated reciprocal classification vignette examples to use `reduction = FALSE`, avoiding unsupported dimensional-reduction usage for `recip_classif` methods.
* Standardized function and parameter naming to snake_case across the active package API.
* Curated naming/style consistency in active source, including replacement of remaining `1:...` indexing patterns.

### Changes

* Removed `get.optimal_clustering` from the public package workflow/documentation while it remains in development.

## Version 0.99.21 (2024-01-27)

### New Features

* Initial Bioconductor submission with ground-truth-agnostic cell type evaluation framework
* Support for multiple input formats: Seurat, SingleCellExperiment, and count matrices
* Multiple dissimilarity methods:
  - Pseudobulk-based distances (Euclidean, Cosine, Pearson)
  - Wasserstein distance on single-cell distributions
  - Reciprocal classification approaches
* Comprehensive internal validation metrics:
  - Silhouette scores
  - Neighborhood Purity
  - Ward's clustering consistency
  - Orbital medoid distances
  - Average similarity measures
  - 2-label silhouette analysis
* Gene selection methods:
  - Highly variable genes (HVG) detection
  - Cell-type-specific marker identification
  - Custom gene list support
* Dimensional reduction support (PCA, pre-computed embeddings)
* Cross-sample and cross-study benchmarking capabilities
* Customizable visualization tools (heatmaps, MDS, PCA plots)
* Hierarchical clustering analysis

### Improvements

* Comprehensive documentation with vignettes
* Examples for all 19 exported functions
* Support for optional dependencies (transformGamPoi, glmGamPoi)
* Efficient euclidean distance computation with custom C implementation
* Parallel processing support via BiocParallel

### Bug Fixes

* Proper handling of sparse matrix formats
* Robust batch effect handling
* Consistent results across different label granularities

### Documentation

* Main vignette: comprehensive tutorial with real-world examples
* Quick start guide: minimal workflow for rapid evaluation
* Full API documentation for all exported functions
