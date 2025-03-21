% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{add.GeneMarkers}
\alias{add.GeneMarkers}
\title{Add Gene Markers to Single-Cell Evaluation Object}
\usage{
add.GeneMarkers(
  scTypeEval,
  ident = NULL,
  sample = NULL,
  method = c("scran.findMarkers"),
  ngenes.total = 500,
  ngenes.celltype = 50,
  ncores = 1,
  bparam = NULL,
  progressbar = TRUE,
  ...
)
}
\arguments{
\item{scTypeEval}{An scTypeEval object.}

\item{ident}{A character string specifying the column name in \code{scTypeEval@metadata} that contains cell type or annotation labels.
If NULL, the function defaults to \code{scTypeEval@active.ident}.}

\item{sample}{Optional. A character string specifying the column name in \code{scTypeEval@metadata} that contains sample identifiers.
If provided, it is used as a blocking factor for differential expression analysis.}

\item{method}{A character string specifying the marker gene identification method. Options are:
\itemize{
\item \code{"scran.findMarkers"} - Uses Scran's \link[scran]{findMarkers} to find differentially expressed genes .
}}

\item{ngenes.total}{Integer specifying the total number of marker genes to retain (default: 500).}

\item{ngenes.celltype}{Integer specifying the number of marker genes per cell type (default: 50).}

\item{ncores}{Integer specifying the number of cores to use for parallel processing (default: 1).}

\item{bparam}{Optional. A BiocParallel parameter object for controlling parallel computation. If provided, overrides \code{ncores}.}

\item{progressbar}{Logical. Whether to display a progress bar (default: TRUE).}

\item{...}{Additional arguments passed to \code{"scran.findMarkers"}.}
}
\value{
The modified \code{scTypeEval} object with gene markers added to \code{scTypeEval@gene.lists[[method]]}.
}
\description{
This function identifies and assigns gene markers for cell types in a single-cell dataset
using \code{scran.findMarkers} method.
}
\details{
\itemize{
\item If \code{ident} is NULL, the function uses \code{scTypeEval@active.ident} as the grouping factor.
\item If \code{sample} is provided, it is used as a blocking factor to adjust for batch effects.
\item The function selects marker genes using \code{scran.findMarkers}.
\item The identified markers are stored in \code{scTypeEval@gene.lists[[method]]}.
}
}
\examples{
\dontrun{
  # Example usage with scran.findMarkers
  sceval <- add.GeneMarkers(sceval, ident = "cell_type", method = "scran.findMarkers")
}

}
