# scTypeEval Unit Tests

This directory contains comprehensive Bioconductor-compliant unit tests for the scTypeEval package.

## Structure

```
tests/
├── testthat.R                               # Test runner
├── README.md                                # This file
├── TEST_SUMMARY.md                          # Summary of all tests
└── testthat/
    ├── helper-data.R                        # Helper functions for test data generation
    ├── test-create-scTypeEval.R             # Tests for scTypeEval object creation
    ├── test-Run-ProcessingData.R            # Tests for data processing
    ├── test-Run-HVG.R                       # Tests for highly variable genes
    ├── test-Run-GeneMarkers.R               # Tests for marker gene identification
    ├── test-Run-PCA.R                       # Tests for PCA computation
    ├── test-Run-Dissimilarity.R             # Tests for dissimilarity computation
    ├── test-get-Consistency.R               # Tests for consistency metrics
    ├── test-utility-functions.R             # Tests for utility functions
    ├── test-S4-classes.R                    # Tests for S4 class structure
    ├── test-integration-workflows.R         # Tests for complete analysis workflows
    ├── test-plotting-functions.R            # Tests for visualization functions
    ├── test-get-NN-hierarchy.R              # Tests for KNN and hierarchical clustering
    └── test-load-SingleCell-object.R        # Tests for data loading utilities
```

## Test Coverage

### Summary

The test suite contains **311 comprehensive tests** across 13 test files:

| Test File | Tests | Functions Covered |
|-----------|-------|-------------------|
| test-create-scTypeEval.R | 17 | create.scTypeEval() |
| test-Run-ProcessingData.R | 24 | Run.ProcessingData() |
| test-Run-HVG.R | 16 | Run.HVG() |
| test-Run-GeneMarkers.R | 15 | Run.GeneMarkers() |
| test-Run-PCA.R | 15 | Run.PCA() |
| test-Run-Dissimilarity.R | 22 | Run.Dissimilarity() |
| test-get-Consistency.R | 20 | get.Consistency() |
| test-utility-functions.R | 18 | set.activeIdent(), add.GeneList(), add.DimReduction(), Add.ProcessedData() |
| test-S4-classes.R | 16 | S4 class structure validation |
| test-integration-workflows.R | 15 | Complete analysis pipelines |
| test-plotting-functions.R | 39 | plot.PCA(), plot.MDS(), plot.Heatmap() |
| test-get-NN-hierarchy.R | 40 | get.NN(), get.hierarchy() |
| test-load-SingleCell-object.R | 54 | load_singleCell_object() |

### Test Details

#### 1. Object Creation (`test-create-scTypeEval.R`; 17 tests)
- Matrix input types (sparse, dense, NULL)
- Seurat object input
- SingleCellExperiment object input
- Metadata validation
- Gene lists and black lists
- Active identity setting

#### 2. Data Processing (`test-Run-ProcessingData.R`; 24 tests)
- Single-cell aggregation
- Pseudobulk aggregation
- Normalization methods
- Filtering (min.samples, min.cells)
- Sample and identity validation
- Verbose output control

#### 3. Highly Variable Genes (`test-Run-HVG.R`; 16 tests)
- `scran` method
- `basic` method
- Sample blocking
- Gene filtering
- Black list handling
- Multiple samples support

#### 4. Marker Genes (`test-Run-GeneMarkers.R`; 15 tests)
- `scran.findMarkers` method
- Gene filtering per cell type
- Black list handling
- Multi-sample support

#### 5. Dimensional Reduction (`test-Run-PCA.R`; 15 tests)
- PCA computation
- Custom gene lists
- Black list handling
- Feature loadings
- Multiple aggregation types

#### 6. Dissimilarity Computation (`test-Run-Dissimilarity.R`; 22 tests)
- Pseudobulk methods (Euclidean, Cosine, Pearson)
- Wasserstein distance
- Reciprocal classification (Match, Score)
- Reduction vs. direct computation
- Gene list filtering
- Black list handling

#### 7. Consistency Metrics (`test-get-Consistency.R`; 20 tests)
- Silhouette
- 2-label silhouette
- Neighborhood purity
- Ward proportion match
- Orbital medoid
- Average similarity
- Multiple metrics simultaneously
- Return format options

#### 8. Utility Functions (`test-utility-functions.R`; 18 tests)
- `set.activeIdent()`
- `add.GeneList()`
- `add.DimReduction()`
- `Add.ProcessedData()`
- Function integration

#### 9. S4 Classes (`test-S4-classes.R`; 16 tests)
- scTypeEval class structure
- DataAssay class structure
- DissimilarityAssay class structure
- DimRed class structure
- Slot validation
- Data integrity

#### 10. Integration Workflows (`test-integration-workflows.R`; 15 tests)
- Complete analysis pipelines
- Multi-step processing chains
- Parameter combinations
- Output validation

#### 11. Plotting Functions (`test-plotting-functions.R`; 39 tests)
- `plot.PCA()` with various parameters and customizations
- `plot.MDS()` with various parameters and customizations
- `plot.Heatmap()` with various parameters and customizations
- Color schemes and aesthetic options
- Legend configurations
- Different aggregation methods
- Dimension selection

#### 12. Nearest Neighbor and Hierarchy (`test-get-NN-hierarchy.R`; 40 tests)
- `get.hierarchy()` with various clustering methods (ward.D2, complete, average, single)
- `get.NN()` with varying KNN parameters (k=3, 5, 10)
- Normalization options
- Return format validation
- Integration with other analysis functions
- Cell type composition analysis

#### 13. Data Loading (`test-load-SingleCell-object.R`; 54 tests)
- Loading Seurat objects (.rds) with split/unsplit modes
- Loading SingleCellExperiment objects (.rds) with split/unsplit modes
- Loading AnnData objects (.h5ad) with counts layer
- Loading AnnData objects (.h5ad) with X layer only
- Loading AnnData objects (.h5ad) with custom layers
- Error handling for missing files and unsupported formats
- Error handling for missing packages
- Data integrity and format validation
- Metadata preservation
- Integration with scTypeEval object creation

## Test Data

The `helper-data.R` file provides functions to generate synthetic single-cell RNA-seq data:

### Main Functions

- **`generate_test_data()`**: Creates a full test dataset with:
  - 500 genes (default)
  - 6 samples (default)
  - Multiple cell types (A, B, C, D)
  - Realistic expression patterns
  - Batch effects
  - Sample-level variation

- **`generate_small_test_data()`**: Creates a smaller dataset for quick tests:
  - 200 genes
  - 5 samples
  - 3 cell types
  - Faster execution

- **`create_test_seurat()`**: Creates Seurat object (if available)
- **`create_test_sce()`**: Creates SingleCellExperiment object (if available)
- **`create_test_scTypeEval()`**: Creates basic scTypeEval object
- **`create_processed_scTypeEval()`**: Creates fully processed scTypeEval object

### Data Characteristics

The synthetic data includes:
- Cell type-specific marker genes
- Background expression
- Sample effects (batch effects)
- Multiple metadata columns (celltype, sample, batch, condition)
- Realistic count distributions using Poisson sampling

## Running Tests

### Quick start
- Full suite: `devtools::test()`
- Single file: `devtools::test(filter = "create-scTypeEval")` (replace with your target file stem)
- Faster reruns in dev: `testthat::test_file("tests/testthat/test-load-SingleCell-object.R")`
- H5AD checks: ensure Python `anndata` is installed (`install.packages("anndata")` will pull via reticulate)

### Run all tests
```r
# From R console
devtools::test()

# Or using testthat directly
testthat::test_dir("tests/testthat")
```

### Run specific test file
```r
testthat::test_file("tests/testthat/test-create-scTypeEval.R")
```

### Run with coverage report
```r
covr::package_coverage()
```

## Bioconductor Compliance

These tests follow Bioconductor standards:

1. **Test Organization**: Uses testthat framework (edition 3)
2. **Helper Functions**: Separate helper file for test data generation
3. **Multiple Samples**: All tests use datasets with ≥5 samples
4. **Parameter Coverage**: Tests cover all function parameters
5. **Edge Cases**: Tests include error conditions and edge cases
6. **S4 Classes**: Tests validate S4 class structure and slots
7. **Dependencies**: Optional tests skip gracefully if dependencies unavailable
8. **Documentation**: All test functions are clearly documented

## Dependencies

### Required for all tests:
- testthat (>= 3.0.0)
- Matrix
- scran
- BiocParallel

### Optional (tests skip if unavailable):
- Seurat
- SingleCellExperiment
- SingleR

## Notes for Developers

### Adding New Tests

1. Create a new test file: `test-<function-name>.R`
2. Use helper functions from `helper-data.R` for test data
3. Test all function parameters
4. Include error conditions
5. Use `skip_if_not_installed()` for optional dependencies
6. Follow naming convention: `test_that("function does X", { ... })`

### Test Best Practices

- Use `verbose = FALSE` in function calls to reduce output
- Test with multiple samples (≥5 recommended)
- Test both success and failure cases
- Use descriptive test names
- Group related tests in the same file
- Keep tests independent and reproducible

### Debugging Tests

```r
# Run a single test interactively
testthat::test_file("tests/testthat/test-create-scTypeEval.R", 
                    reporter = "location")

# Run with detailed output
devtools::test(reporter = "check")
```

## Continuous Integration

These tests are designed to run in CI/CD pipelines:
- Compatible with GitHub Actions
- Compatible with Bioconductor build system
- Tests complete in reasonable time
- Graceful handling of missing optional dependencies

## Contact

For questions about the tests, please open an issue at:
https://github.com/carmonalab/scTypeEval/issues
