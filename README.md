# scTypeEval :microscope: :chart_with_upwards_trend:

<p align="center">

<img src="docs/scTypeEval_logo.png" height="140"/>

</p>

## A Framework for Evaluating Single-Cell Transcriptomics Cell Type Annotations

`scTypeEval` is a comprehensive tool to measure the quality of cell type annotations in single-cell RNA-seq data. It provides a suite of validation metrics to benchmark annotations and improve classification performance.

Accurately evaluating cell type annotations in scRNA-seq is challenging due to the lack of actual ground truth. Traditional approaches — including manual curation and marker-based validation — are subjective, inconsistent, and often overlook subtle expression patterns. However, reliable cell type annotation should be consistent across samples and datasets, ensuring reproducibility and biological relevance.

scTypeEval addresses these challenges by using internal validation metrics to assess label consistency without requiring predefined annotations. Tested across diverse datasets and classification scenarios, it provides a robust, unbiased framework for benchmarking annotation performance.

### Key features of scTypeEval

1. **Internal validation metrics**: Assess the quality of cell type annotations without requiring ground truth labels.
2. **Benchmarking across datasets**: Evaluate annotation performance across different samples and any classification scenarios (manual annotation, automatic classifiers or clustering).
3. **Customizable framework**: Support for user-defined different gene lists, normalization methods, and metrics parameters.
4. **Robust evaluation**: Identify misclassifications and compare classification methods using diverse metrics.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("carmonalab/scTypeEval")
```


## Usage



### Main consistency metrics
Summary of the main consistency metrics available in `scTypeEval`:

| Metric              | Description                                                            | Type                           | Quantification                                                        | Level      | Hyperparameters                |
|---------------------|------------------------------------------------------------------------|-------------------------------|------------------------------------------------------------------------|------------|--------------------------------|
| Silhouette         | Contrast intra-cluster tightness with inter-cluster separation        | Distance-based                | Mean silhouette width per cell type                                   | Cell type  |                                |
| Modularity         | Strength of intra-community density compared to random expectation   | Graph structure               | Global modularity score and per-cell-type modularity contribution    | Global     | Number of neighbors            |
| Ward              | Hierarchical clustering minimizing intra-cluster variance             | Clustering-based              | Proportion Match: Normalized proportion of reference labels in the dominant cluster | Cell type  | Linkage method (e.g., ward.D2) |
|                   |                                                                        |                               | NMI: Normalized shared information between cluster assignments and reference labels | Global     |                                |
| Leiden            | Community detection for well-connected communities                    | Clustering-based              | ARI: Adjusted overlap of cluster assignments with reference labels   | Global     | Number of neighbors            |
| Neighborhood Purity | Proportion of K-Nearest Neighbors sharing the same label            | Local agreement               | Mean proportion of neighbors sharing the same identity, normalized by expected proportions | Cell type  | Number of neighbors           |
| Mutual Best Hit   | Consistency of annotations by applying a classifier in both directions | Inter-sample similarity      | Score: product of reciprocal prediction scores                        | Cell type  | Classifier (e.g. SingleR)      |
|                   |                                                                        |                               | Best Match: Normalized proportion of reciprocal prediction            |            |                                |
| Graph Connectivity | Proportion of elements within a group that belong to its largest connected component | Graph structure | Size of the largest connected component relative to reference group size | Cell type  | Number of neighbors            |



## scTypeEval demos

Find a vignette describing scTypeEval main functions in [html]() and its [code (repository)](https://github.com/carmonalab/scTypeEval_CaseStudies).
