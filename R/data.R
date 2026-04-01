#' Default Gene Blacklist for scTypeEval
#'
#' @name black_list
#'
#' @description
#' A curated list of genes typically excluded from single-cell RNA-seq analysis
#' to reduce technical artifacts and improve cell type annotation quality.
#'
#' @format A character vector with 7,775 gene symbols including:
#' \describe{
#'   \item{Mitochondrial genes}{MT- prefix genes (MT-RNR1, MT-RNR2, MT-TA, etc.)}
#'   \item{Ribosomal genes}{RPL and RPS prefix genes}
#'   \item{Long non-coding RNAs}{Genes with -AS1, -AS2, -IT1, -DT suffixes}
#'   \item{MicroRNAs}{MIR prefix genes (MIR1-1, MIR1-2, etc.)}
#'   \item{Small nuclear RNAs}{SNORA, SNORD, RNU prefix genes}
#'   \item{Cell cycle genes}{G1/S and G2/M phase markers}
#'   \item{Sex chromosome genes}{X and Y chromosome specific genes}
#'   \item{Other technical artifacts}{Heat shock proteins, immediate early genes}
#' }
#'
#' @details
#' This blacklist is automatically loaded and used by default in scTypeEval
#' when \code{black_list = NULL} is specified in functions like \code{\link{run_hvg}}
#' and \code{\link{run_gene_markers}}. Users can override this default by providing
#' a custom character vector of gene symbols to exclude.
#'
#' The blacklist helps improve downstream analysis by removing:
#' \itemize{
#'   \item Genes with high technical variance unrelated to cell type identity
#'   \item Genes that may confound clustering (e.g., cell cycle genes)
#'   \item Genes with batch-specific expression patterns
#'   \item Non-coding RNAs that may not be informative for cell type annotation
#' }
#'
#' @return A character vector containing gene symbols to be excluded from analysis,
#' including cell cycle genes (G1/S and G2/M), mitochondrial genes, ribosomal genes,
#' TCR and immunoglobulin genes, pseudogenes, heat shock proteins, non-coding RNAs,
#' and sex chromosome genes (X and Y).
#'
#' @source
#' Generated using the SignatuR package (Andreatta et al., \url{https://github.com/carmonalab/SignatuR}) with the following categories:
#' \itemize{
#'   \item Cell cycle genes: G1/S and G2/M phase markers from SignatuR::Programs
#'   \item Technical artifacts: Mitochondrial, ribosomal, TCR, immunoglobulins, pseudogenes, HSP, and non-coding RNAs from SignatuR::Blocklists
#'   \item Sex chromosome genes: X-inactivation escapees and Y-chromosome specific genes from GenderGenes database (\url{http://bioinf.wehi.edu.au/software/GenderGenes/})
#' }
#'
#' @examples
#' # Load the default gene blacklist
#' data(black_list)
#'
#' # Inspect the first few entries
#' head(black_list)
#'
#' @seealso \code{\link{run_hvg}}, \code{\link{run_gene_markers}}
#'
#' @keywords datasets
NULL
