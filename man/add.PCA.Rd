% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{add.PCA}
\alias{add.PCA}
\title{Perform PCA on a Gene List and Store Results in scTypeEval object}
\usage{
add.PCA(
  scTypeEval,
  ident = NULL,
  sample = NULL,
  normalization.method = c("Log1p", "CLR", "pearson"),
  gene.list = NULL,
  data.type = c("sc", "pseudobulk", "pseudobulk_1vsall"),
  min.samples = 5,
  min.cells = 10,
  black.list = NULL,
  ndim = 30,
  ncores = 1,
  bparam = NULL,
  progressbar = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{scTypeEval}{A \code{scTypeEval} object containing single-cell expression data.}

\item{ident}{Character. Metadata column name used to group cells (e.g., cell type annotation).
If \code{NULL}, the active identity from \code{scTypeEval} is used.}

\item{sample}{Character. Metadata column name specifying sample identity for pseudobulk analysis.
Required for pseudobulk and \code{pseudobulk_1vsall} data types.}

\item{normalization.method}{Character. Method for normalizing gene expression before PCA. See \link[scTypeEval]{add.HVG} for more details.
Options: \code{"Log1p"}, \code{"CLR"}, \code{"pearson"} (default: \code{"Log1p"}).}

\item{gene.list}{Named list of character vectors. Each element is a set of genes for PCA analysis.
If \code{NULL}, all pre-defined gene lists in \code{scTypeEval} are used recursively.}

\item{data.type}{Character. Type of data to analyze. Options: \code{"sc"}, \code{"pseudobulk"},
or \code{"pseudobulk_1vsall"}. Default is \code{"sc"}.}

\item{min.samples}{Integer. Minimum number of samples required for pseudobulk PCA (default: 5).}

\item{min.cells}{Integer. Minimum number of cells required per group for PCA (default: 10).}

\item{black.list}{Character vector. Genes to exclude from PCA. If \code{NULL}, defaults to
the blacklist stored in \code{scTypeEval}.}

\item{ndim}{Integer. Number of principal components to compute (default: 30).}

\item{ncores}{Integer. Number of CPU cores to use for parallel processing. Default: \code{1}.}

\item{bparam}{Parallel backend parameter object for BiocParallel. If provided, overrides \code{ncores}.}

\item{progressbar}{Logical. Whether to display a progress bar during computation. Default: \code{FALSE}.}

\item{verbose}{Logical. Whether to print messages during execution. Default: \code{TRUE}.}
}
\value{
The modified \code{scTypeEval} object with PCA results stored in \code{reductions} slot.
}
\description{
This function computes Principal Component Analysis (PCA) based on an specified gene list
and stores the results in the \code{reductions} slot of an scTypeEval object.
}
\examples{
\dontrun{
sceval <- add.PCA(sceval, # scTypeEval object
                 ident = "cell_type",
                 sample = "sample_id",
                 data.type = "pseudobulk")
}

}
\seealso{
\link{add.HVG}
}
