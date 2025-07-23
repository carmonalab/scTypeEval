# Create scTypeEval object

#' @title Create an scTypeEval object for evaluating cell type classification.
#'
#' @description This function initializes an `scTypeEval` object from various input formats, including Seurat, 
#' SingleCellExperiment, or raw count matrices. It ensures compatibility by validating input types 
#' and structures before constructing the object.
#'
#' @param matrix A Seurat object, SingleCellExperiment object, or a count matrix (dense or sparse).
#' @param metadata A metadata dataframe (required if `matrix` is a raw matrix). 
#'        It must have the same number of rows as columns in the count matrix, sharing rownames and colnames respectively.
#' @param gene.lists A named list of gene sets to use in the evaluation.
#' @param black.list A character vector of genes to exclude from analysis.
#' @param active.ident The active identity class or cluster label (optional).
#' @param version A string indicating the version of the `scTypeEval` object (optional).
#'
#' @return An `scTypeEval` object containing:
#' \itemize{
#'   \item \code{counts}: A sparse count matrix (dgCMatrix).
#'   \item \code{metadata}: A dataframe with metadata for each cell.
#'   \item \code{gene.lists}: A list of gene sets used in classification.
#'   \item \code{black.list}: A vector of excluded genes.
#'   \item \code{active.ident}: Active cluster identity (if provided).
#'   \item \code{version}: The version tag for the object.
#' }
#'
#' @details The function converts input objects into a standardized format before constructing 
#'          an `scTypeEval` object. It ensures that the count matrix and metadata are correctly 
#'          aligned and validates input types before proceeding.
#'
#' @examples
#' \dontrun{
#' # From count matrix and metadata dataframe
#' sceval <- create.scTypeEval(matrix = count_matrix, metadata = metadata)
#' 
#' # From Seurat object
#' sceval <- create.scTypeEval(seurat_obj)
#' }
#' 
#' @importClassesFrom Matrix dgCMatrix
#'
#' @export create.scTypeEval


create.scTypeEval <- function(matrix = NULL,
                              metadata = NULL,
                              gene.lists = list(),
                              black.list = character(),
                              active.ident = NULL,
                              version = "2.0.0") {
   
   # Handle NULL matrix case
   if (is.null(matrix)) {
      message("⚠️  No matrix provided. Initializing scTypeEval with empty counts.")
      
      # Create an empty dgCMatrix
      counts <- Matrix::Matrix(0, nrow = 0, ncol = 0, sparse = TRUE)
      
      # If no metadata provided, create an empty data frame
      if (is.null(metadata)) {
         stop("metadata dataframe must be provided.\n")
      } else {
         metadata <- as.data.frame(metadata)
      }
      
   } else if (inherits(matrix, "Seurat")) {
      counts <- as(matrix@assays$RNA@counts, "dgCMatrix")
      metadata <- as.data.frame(matrix@meta.data)
      
   } else if (inherits(matrix, "SingleCellExperiment")) {
      counts <- as(SummarizedExperiment::assay(matrix, "counts"), "dgCMatrix")
      metadata <- as.data.frame(SummarizedExperiment::colData(matrix))
      
   } else if (inherits(matrix, "matrix") || inherits(matrix, "dgCMatrix")) {
      counts <- as(matrix, "dgCMatrix")
      if (is.null(metadata)) {
         stop("For matrix input, metadata dataframe must be provided.")
      }
      metadata <- as.data.frame(metadata)
      
   } else {
      stop("Input object must be NULL, Seurat, SingleCellExperiment, or matrix-like object.")
   }
   
   # Check dimensions (skip if matrix is NULL/empty)
   if (!is.null(matrix) && ncol(counts) != nrow(metadata)) {
      stop("Different number of columns in counts and number of rows in metadata.")
   }
   
   # Create the scTypeEval object
   scTypeEval_obj <- methods::new("scTypeEval",
                                  counts = counts,
                                  metadata = metadata,
                                  gene.lists = gene.lists,
                                  black.list = black.list,
                                  active.ident = active.ident,
                                  version = version)
   
   tmp <- .check_ident(scTypeEval_obj,
                       active.ident,
                       verbose = FALSE)
   
   return(scTypeEval_obj)
}


#' @title Set the active ident or annotation labels for an scTypeEval object.
#'
#' @description This function assigns a specific cell type annotation from the metadata as the active identity 
#' in an `scTypeEval` object, allowing for downstream analysis based on the selected classification.
#'
#' @param scTypeEval An `scTypeEval` object.
#' @param ident A character string specifying the column in `scTypeEval@metadata` to set as the active identity.
#'
#' @return The modified `scTypeEval` object with the active identity set.
#'
#' @details The function ensures that the provided identity exists within the metadata before setting it. 
#'          If no valid identity is provided, an error is raised.
#'
#' @examples
#' \dontrun{
#' sceval <- create.scTypeEval(count_matrix, metadata = metadata)
#' sceval <- set.activeIdent(sceval, ident = "cell_type")
#' }
#'
#' @export set.activeIdent


set.activeIdent <- function(scTypeEval,
                            ident = NULL){
   if(is.null(ident)){
      stop("Specificy a cell type annotation in the provided metadata")
   }
   if(!ident %in% names(scTypeEval@metadata)){
      stop("Please provide a ident, i.e. a cell type or annotation
           to group cells included in metadata")
   }
   
   scTypeEval@active.ident <- ident
   
   return(scTypeEval)
}


Run.ProcessingData <- function(scTypeEval,
                               ident = NULL,
                               sample = NULL,
                               normalization.method = "Log1p",
                               min.samples = 5,
                               min.cells = 10,
                               ndim = 30,
                               verbose = TRUE){
   #### preflights checks
   ident.name <- ident
   ident <- .check_ident(scTypeEval, ident, verbose = verbose)
   sample <- .check_sample(scTypeEval, sample, verbose = verbose)
   
   # set normalization method
   normalization.method <- normalization.method[1]
   
   #### RUN 
   data.list <- lapply(aggregation_types,
                       function(ag){
                          if(verbose){message("# Processing data for ", ag, " ... \n")}
                          # Transform data
                          if(verbose){message("   Transforming and filtering count matrix... \n")}
                          mat <- get.matrix(scTypeEval@counts,
                                            ident = ident,
                                            sample = sample,
                                            aggregation = aggregation,
                                            min.samples = min.samples,
                                            min.cells = min.cells)
                          
                          # Normalize data
                          if(verbose){message("   Normalizing count matrix... \n")}
                          norm.mat <- Normalize_data(mat@matrix,
                                                     method = normalization.method)
                          
                          rr <- methods::new("DataAssay",
                                             data = norm.mat,
                                             aggregation = ag,
                                             group = mat@group,
                                             ident = list(ident.name = mat@ident),
                                             sample = mat@sample)
                          return(rr)
                       })
   names(data.list) <- aggregation_types
   
   # add to scTypeEval object
   scTypeEval@data <- data.list
   
   return(scTypeEval)
   
}


Add.ProcessedData <- function(scTypeEval,
                              data,
                              aggregation,
                              ident = NULL,
                              ident.name = "custom",
                              sample = NULL,
                              filter = FALSE,
                              min.samples = 5,
                              min.cells = 10,
                              verbose = TRUE){
   
   # preflight
   if(!aggregation %in% aggregation_types){
      stop("Only supported aggregations are just `single-cell` (no aggregation) or `pseudobulk` per cell type and sample.\n")
   }
   #### preflights checks
   ident <- .check_ident(ident = ident, verbose = verbose)
   if(length(ident) == ncol(data)){
      stop("Different length of ident and ncol of data\n")
   }
   sample <- .check_sample(sample = sample, verbose = verbose)
   if(length(sample) == ncol(data)){
      stop("Different length of sample and ncol of data\n")
   }
   groups <- factor(paste(sample, ident, sep = "_"))
   
   if(aggregation == "single-cell"){
      
      if(filter){
         if(verbose){message(" - Filtering cells based on min.samples and min.cells\n")}
         
         mat <- get.matrix(data,
                           ident = ident,
                           sample = sample,
                           aggregation = aggregation,
                           min.samples = min.samples,
                           min.cells = min.cells)
         
         data <- mat@matrix
         ident <- mat@ident
         sample <- mat@sample
         groups <- factor(paste(sample, ident, sep = "_"))
         
      } 
      
   } else if(aggregation == "pseudobulk"){
      if(filter) {
         stop("Filtering only allowed with single-cell data aggregation\n.")
      }
      if(verbose){message("Expecting cells aggregated (pseudobulk) per cell type (ident) and sample of origin.\n")}
      
      prop <- table(groups)
      if(any(prop)>1){
         stop("Aggregated data by pseudobulk expect a frequency of 1 for each cell type & sample combination.\n")
      }
   }
   
   
   # Create Data assay
   rr <- methods::new("DataAssay",
                      data = data,
                      aggregation = aggregation,
                      group = group,
                      ident = list(ident.name = ident),
                      sample = sample)
   
   # add to scTypeEval object
   scTypeEval@data[[aggregation]]<- rr
   
   return(scTypeEval)
   
}

#' @title Identify and add highly variable genes (HVG) to an scTypeEval object.
#'
#' @description This function detects highly variable genes in an `scTypeEval` object using different 
#' normalization and variance calculation methods. The identified genes are stored in the 
#' `gene.lists` slot under "HVG".
#'
#' @param scTypeEval An `scTypeEval` object.
#' @param normalization.method Character string specifying the normalization method to apply.
#'                              Options: `"Log1p"`, `"CLR"`, `"pearson"` (default: `"Log1p"`).
#' @param var.method Character string specifying the method for identifying highly variable genes.
#'                   Options: `"scran"` or `"basic"` (default: `"scran"`).
#' @param sample Optional. A metadata column name to use for grouping cells before computing HVGs.
#' @param ngenes Integer specifying the number of highly variable genes to retain (default: `500`).
#' @param black.list A character vector of genes to exclude from HVG selection (default: `NULL`, uses `scTypeEval@black.list`).
#' @param ncores Integer specifying the number of CPU cores to use for parallel processing (default: `1`).
#' @param bparam A `BiocParallel` backend parameter object for parallelization. If provided, overrides `ncores`.
#' @param progressbar Logical, whether to display a progress bar during computation (default: `TRUE`).
#' @param verbose Logical. Whether to print messages during execution. Default: `TRUE`.
#' @param ... Additional arguments passed to pearson residuals normalization functions.
#'
#' @return The modified `scTypeEval` object with HVGs added to `scTypeEval@gene.lists[["HVG"]]`.
#'
#' @details The function first normalizes the expression matrix using the specified method:  
#' - `"Log1p"`: Applies log(x+1) transformation to stabilize variance in the count matrix. Recommended for scran `modelGeneVar`.
#' - `"CLR"`: Performs centered log-ratio normalization, where each gene’s expression  
#'            is divided by the geometric mean across all genes in a given cell.  
#' - `"pearson"`: Computes Pearson residuals from a negative binomial model  
#'            fitted to the count data. This method is commonly used to reduce technical noise  
#'            while preserving biological signal. Not recommended to get HVG, as variance stabilization may suppress variability.
#'
#' After normalization, HVGs are computed using one of the following approaches:  
#' - `"scran"`: Uses the **scran** package’s \link[scran]{modelGeneVar} function to identify highly variable genes (HVGs).  
#'            This approach applies model-based variance decomposition to detect genes with significant biological variation  
#'            while accounting for technical noise, sample, and batch effects.  
#'  - `"basic"`: A simpler method based on the variance-to-mean ratio. It ranks genes  
#'            by their coefficient of variation (CV) and selects the most variable ones.  
#'            This approach is computationally efficient but lacks statistical modeling.  
#' Any genes in the `black.list` is removed before HVG selection.
#'
#' @examples
#' \dontrun{
#' # create scTypeEval object
#' sceval <- create.scTypeEval(count_matrix, metadata = metadata)
#' # Compute and add HVG
#' sceval <- add.HVG(sceval)
#' }
#'
#' @export Run.HVG




Run.HVG <- function(scTypeEval,
                    normalization.method = c("Log1p", "CLR", "pearson"),
                    var.method = c("scran", "basic"),
                    sample = NULL,
                    ngenes = 2000,
                    black.list = NULL,
                    ncores = 1,
                    bparam = NULL,
                    progressbar = FALSE,
                    verbose = TRUE,
                    ...){
   
   black.list <- .check_blacklist(scTypeEval, black.list, verbose = verbose)
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   
   # normalized matrix
   norm.mat <- scTypeEval@data[["single-cell"]]
   if(is.null(norm.mat)){
      stop("No normalization slot found. Please run before `Run.ProcessingData()`.\n")
   }
   mat <- norm.mat@matrix
   if(!is.null(sample)){
      sample <- norm.mat@sample
   }
   
   # remove blacked listed genes
   if(!is.null(black.list) && verbose){message("Filtering out black listed genes... \n")}
   mat <- mat[!rownames(mat) %in% black.list,]
   
   var.method <- var.method[1]
   
   # get highly variable genes
   if(verbose){message("Computing HVG... \n")}
   hgv <- switch(var.method,
                 "basic" = get.HVG(norm.mat = mat,
                                   ngenes = ngenes,
                                   sample = sample,
                                   bparam = param),
                 "scran" = get.GeneVar(norm.mat = mat,
                                       sample = sample,
                                       ngenes = ngenes,
                                       bparam = param),
                 stop(var.method, " not supported for getting variable genes.")
   )
   
   
   scTypeEval@gene.lists[["HVG"]] <- hgv
   
   return(scTypeEval)
}


#' @title Add Gene Markers to Single-Cell Evaluation Object
#'
#' @description This function identifies and assigns gene markers for cell types in a single-cell dataset 
#' using `scran.findMarkers` method.
#'
#' @param scTypeEval An scTypeEval object.
#' @param ident A character string specifying the column name in `scTypeEval@metadata` that contains cell type or annotation labels. 
#'   If NULL, the function defaults to `scTypeEval@active.ident`.
#' @param sample Optional. A character string specifying the column name in `scTypeEval@metadata` that contains sample identifiers.
#'   If provided, it is used as a blocking factor for differential expression analysis.
#' @param method A character string specifying the marker gene identification method. Options are:
#'   \itemize{
#'     \item `"scran.findMarkers"` - Uses Scran's \link[scran]{findMarkers} to find differentially expressed genes .
#'   }
#' @param ngenes.total Integer specifying the total number of marker genes to retain (default: 500).
#' @param ngenes.celltype Integer specifying the number of marker genes per cell type (default: 50).
#' @param ncores Integer specifying the number of cores to use for parallel processing (default: 1).
#' @param bparam Optional. A BiocParallel parameter object for controlling parallel computation. If provided, overrides `ncores`.
#' @param progressbar Logical. Whether to display a progress bar (default: TRUE).
#' @param verbose Logical. Whether to print messages during execution. Default: `TRUE`.
#' @param ... Additional arguments passed to `"scran.findMarkers"`.
#'
#' @return The modified `scTypeEval` object with gene markers added to `scTypeEval@gene.lists[[method]]`.
#' 
#' @details 
#' - If `ident` is NULL, the function uses `scTypeEval@active.ident` as the grouping factor.
#' - If `sample` is provided, it is used as a blocking factor to adjust for batch effects.
#' - The function selects marker genes using `scran.findMarkers`.
#' - The identified markers are stored in `scTypeEval@gene.lists[[method]]`.
#'
#' @examples
#' \dontrun{
#'   # Example usage with scran.findMarkers
#'   sceval <- add.GeneMarkers(sceval, ident = "cell_type", method = "scran.findMarkers")
#' }
#'
#' @export add.GeneMarkers


Run.GeneMarkers <- function(scTypeEval,
                            method = c("scran.findMarkers"),
                            ngenes.celltype = 50,
                            black.list = NULL,
                            ncores = 1,
                            bparam = NULL,
                            progressbar = FALSE,
                            verbose = TRUE,
                            ...){
   
   method <- method[1]
   if(!method %in% c("scran.findMarkers")){
      stop("Supported gene markes definitions is `scran.findMarkers`")
   }
   
   black.list <- .check_blacklist(scTypeEval, black.list, verbose = verbose)
   
   # normalized matrix
   norm.mat <- scTypeEval@data[["single-cell"]]
   if(is.null(norm.mat)){
      stop("No normalization slot found. Please run before `Run.ProcessingData()`.\n")
   }
   mat <- norm.mat@matrix
   ident <- unlist(norm.mat@ident)
   sample <- norm.mat@sample
   
   # remove blacked listed genes
   if(!is.null(black.list) && verbose){message("Filtering out black listed genes... \n")}
   norm.mat <- norm.mat[!rownames(norm.mat) %in% black.list,]
   
   if(verbose){message("Computing cell type markers for", names(norm.mat@ident),  "... \n")}
   markers <- switch(method,
                     "scran.findMarkers" = get.DEG(mat = mat,
                                                   ident = ident,
                                                   block = sample,
                                                   ngenes.celltype = ngenes.celltype,
                                                   ncores = ncores,
                                                   bparam = bparam,
                                                   progressbar = progressbar,
                                                   ...)
   )
   
   scTypeEval@gene.lists[[method]] <- markers
   
   
   return(scTypeEval)
   
}

#' @title Add a gene list to an scTypeEval object.
#'
#' @description This function appends a new gene list to an existing `scTypeEval` object.
#' It ensures that the input list is valid and assigns names if they are missing.
#'
#' @param scTypeEval An `scTypeEval` object to which the gene list will be added.
#' @param gene.list A named list of gene sets to append. If the list is unnamed, names will be assigned automatically.
#'
#' @return An updated `scTypeEval` object with the new gene list added.
#'
#' @details The function verifies that `gene.list` is provided and is a valid list.
#'          If any elements in `gene.list` lack names, they are automatically renamed.
#'
#' @examples
#' \dontrun{
#' sceval <- add.GeneList(sceval, gene.list = list("cytokines" = c("IL10", "IL6", "IL4",...)))
#' }
#'
#' @export add.GeneList



add.GeneList <- function(scTypeEval,
                         gene.list = NULL){
   
   # check if it is a list and if it is named
   # Check if gene.list is provided
   if (is.null(gene.list)) {
      stop("gene.list cannot be NULL. Please provide a valid list.")
   }
   
   if(!is.list(gene.list)){
      gene.list <- list(gene.list)
   }
   
   # If gene.list is unnamed, assign names
   if (any(is.null(names(gene.list)))) {
      names(gene.list) <- paste0("gene.list", seq_along(gene.list))
      warning("All or some names of the list is NULL, renaming list.")
   }
   
   scTypeEval@gene.lists <- c(scTypeEval@gene.lists, gene.list)
   
   return(scTypeEval)
}





#' @title Perform PCA on a Gene List and Store Results in scTypeEval object
#'
#' @description This function computes Principal Component Analysis (PCA) based on an specified gene list
#' and stores the results in the \code{reductions} slot of an scTypeEval object.
#'
#' @param scTypeEval A \code{scTypeEval} object containing single-cell expression data.
#' @param ident Character. Metadata column name used to group cells (e.g., cell type annotation).
#' If \code{NULL}, the active identity from \code{scTypeEval} is used.
#' @param sample Character. Metadata column name specifying sample identity for pseudobulk analysis.
#' Required for \code{pseudobulk}, \code{pseudobulk_1vsall}, and \code{GloScope} data types.
#' @param normalization.method Character. Method for normalizing gene expression before PCA. See \link[scTypeEval]{add.HVG} for more details.
#' Options: \code{"Log1p"}, \code{"CLR"}, \code{"pearson"} (default: \code{"Log1p"}).
#' @param gene.list Named list of character vectors. Each element is a set of genes for PCA analysis.
#' If \code{NULL}, all pre-defined gene lists in \code{scTypeEval} are used recursively.
#' @param data.type Character. Type of data to analyze. Options: \code{"sc"}, \code{"pseudobulk"}, or
#' \code{"pseudobulk_1vsall"}. Default is \code{"sc"}.
#' @param min.samples Integer. Minimum number of samples required for pseudobulk PCA (default: 5).
#' @param min.cells Integer. Minimum number of cells required per group for PCA (default: 10).
#' @param black.list Character vector. Genes to exclude from PCA. If \code{NULL}, defaults to
#' the blacklist stored in \code{scTypeEval}.
#' @param ndim Integer. Number of principal components to compute (default: 30).
#' @param ncores Integer. Number of CPU cores to use for parallel processing. Default: `1`.
#' @param bparam Parallel backend parameter object for BiocParallel. If provided, overrides `ncores`.
#' @param progressbar Logical. Whether to display a progress bar during computation. Default: `FALSE`.
#' @param verbose Logical. Whether to print messages during execution. Default: `TRUE`.
#'
#' @return The modified \code{scTypeEval} object with PCA results stored in \code{reductions} slot.

#'
#' @examples
#' \dontrun{
#' sceval <- add.PCA(sceval, # scTypeEval object
#'                  ident = "cell_type",
#'                  sample = "sample_id",
#'                  data.type = "pseudobulk")
#' }
#' 
#' @seealso \link{add.HVG}
#'
#' @export Run.PCA


# obtain PCA from a gene list
Run.PCA <- function(scTypeEval,
                    gene.list,
                    black.list,
                    ndim = 30,
                    verbose = TRUE){
   
   if(length(scTypeEval@data)<1){
      stop("No normalization slot found. Please run before `Run.ProcessingData()`.\n")
   }
   
   gene.list <- .check_genelist(scTypeEval, gene.list, verbose = verbose)
   black.list <- .check_blacklist(scTypeEval, black.list, verbose = verbose)
   
   pca.list <- lapply(names(scTypeEval@data),
                      function(ag){
                         if(verbose){message("# Computing PCA data for ", ag, " ... \n")}
                         # normalized matrix
                         mat <- scTypeEval@data[[ag]]
                         if(is.null(mat)){
                            stop("No normalization slot found. Please run before `Run.ProcessingData()`.\n")
                         }
                         
                         mat <- .general_filtering(mat,
                                                   black.list = black.list,
                                                   gene.list = gene.list,
                                                   verbose = verbose)
                         
                         # compute PCA
                         if(verbose){message("   Computing PCA space... \n")}
                         pr <- custom_prcomp(mat@matrix,
                                             ndim = ndim,
                                             verbose = verbose)
                         
                         # Create DimRed assay
                         rr <- methods::new("DimRed",
                                            embeddings = t(pr$x),
                                            feature.loadings = pr$rotation,
                                            gene.list = gene.list,
                                            black.list = black.list,
                                            aggregation = ag,
                                            group = mat@group,
                                            ident = mat@ident,
                                            sample = mat@sample,
                                            key = "PCA")
                      })
   
   names(pca.list) <- names(scTypeEval@data)
   
   # add to scTypeEval object
   for(n in names(pca.list)){
      scTypeEval@consistency[[n]] <- pca.list[[n]]
   }
   
   return(scTypeEval)
   
}

add.DimReduction <- function(scTypeEval,
                             embeddings,
                             aggregation,
                             ident = NULL,
                             ident.name = "custom",
                             sample = NULL,
                             key = NULL,
                             gene.list = NULL,
                             black.list = NULL,
                             feature.loadings = NULL,
                             filter = FALSE,
                             min.samples = 5,
                             min.cells = 10,
                             verbose = TRUE
                             
)
{
   # preflight
   if(!aggregation %in% aggregation_types){
      stop("Only supported aggregations are just `single-cell` (no aggregation) or `pseudobulk` per cell type and sample.\n")
   }
   ident <- .check_ident(ident = ident, verbose = verbose)
   if(length(ident) == ncol(data)){
      stop("Different length of ident vector and ncol of data\n")
   }
   sample <- .check_sample(sample = sample, verbose = verbose)
   if(length(sample) == ncol(data)){
      stop("Different length of sample vector and ncol of data\n")
   }
   groups <- factor(paste(sample, ident, sep = "_"))
   
   if(aggregation == "single-cell"){
      
      if(filter){
         if(verbose){message(" - Filtering cells based on min.samples and min.cells\n")}
         
         mat <- get.matrix(data,
                           ident = ident,
                           sample = sample,
                           aggregation = aggregation,
                           min.samples = min.samples,
                           min.cells = min.cells)
         
         embeddings <- mat@matrix
         ident <- mat@ident
         sample <- mat@sample
         groups <- factor(paste(sample, ident, sep = "_"))
         
      } 
      
   } else if(aggregation == "pseudobulk"){
      if(filter) {
         stop("Filtering only allowed with single-cell data aggregation\n.")
      }
      if(verbose){message("Expecting cells aggregated (pseudobulk) per cell type (ident) and sample of origin.\n")}
      
      prop <- table(groups)
      if(any(prop)>1){
         stop("Aggregated data by pseudobulk expect a frequency of 1 for each cell type & sample combination.\n")
      }
   }
   
   # Create DimRed assay
   rr <- methods::new("DimRed",
                      embeddings = emdeddings,
                      feature.loadings = feature.loadings,
                      gene.list = gene.list,
                      black.list = black.list,
                      aggregation = aggregation,
                      group = groups,
                      ident = list(ident.name = ident),
                      sample = sample,
                      key = key)
   
   # add to scTypeEval object
   scTypeEval@reductions[[aggregation]]<- rr
   
   return(scTypeEval)
   
}


Run.Dissimilarity <- function(scTypeEval,
                              method = "WasserStein",
                              reduction = TRUE,
                              gene.list = NULL,
                              black.list = NULL,
                              BestHit.classifier = "SingleR",
                              ncores = 1,
                              bparam = NULL,
                              progressbar = FALSE,
                              verbose = TRUE
                              
){
   # run preflight checks
   if(!method %in% names(dissimilarity_methods)){
      stop("Not supported dissimilarity method.\n 
           Please choice one of: ",
           paste(names(dissimilarity_methods), collapse = ", "),
           ".\n")
   }
   
   slot <- dissimilarity_methods[method]
   
   if(reduction && !method %in% no_dr_ds){
      warning("No dimensional reduction dissimilarity computationg supported for: ",
              paste(no_dr_ds, collapse = ", "), "swithcing `reduction=FALSE`")
      reduction <- FALSE
   }
   
   if(reduction){
      mat_ident <- scTypeEval@reductions[[slot]]
      if(is.null(mat_ident)){
         stop("No dimensional reduction slot found for ", slot ,
              ". Please run before `Run.PCA()` or add a valid dimensional reduction assay.\n")
      }
      mat <- mat_ident@embeddings
      gene.list <- mat_ident@gene.list
      black.list <- mat_ident@black.list
      
   } else {
      mat_ident <- scTypeEval@data[[slot]]
      if(is.null(mat_ident)){
         stop("No processed data slot found for ", slot ,
              ". Please run before `Run.ProcessingData()` or add a data assay.\n")
      }
      gene.list <- .check_genelist(scTypeEval, gene.list, verbose = verbose)
      black.list <- .check_blacklist(scTypeEval, black.list, verbose = verbose)
      mat_ident <- .general_filtering(mat_ident,
                                      black.list = black.list,
                                      gene.list = gene.list,
                                      verbose = verbose)
      mat <- mat_ident@data
   }
   
   
   # set paralelization
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   # split matrix per sample for BestHit, and for celltype & sample for Wasserstein, no need for the rest
   
   aggregation <- strsplit(method, ":")[[1]][1] |> tolower()
   dist.type <- strsplit(method, ":")[[1]][2] |> tolower()
   
   dist <- switch(
      aggregation,
      "pseudobulk" = get.distance(norm.mat = mat, distance.method = dist.type, verbose = verbose),
      "wasserstein" = compute_wasserstein(mat = mat, group = mat_ident@group, bparam = param, verbose = verbose),
      "besthit" = bestHit(mat = mat, classifier = BestHit.classifier, method = dist.type, bparam = param, verbose = verbose),
      stop(aggregation, "is not a supported method.\n")
      )
   
   # build dissimilarity object
   rr <- methods::new("DissimilarityAssay",
                      dissimilarity = dist,
                      method = method,
                      gene.list = gene.list,
                      black.list = black.list,
                      aggregation = aggregation,
                      ident = mat_ident@ident,
                      sample = mat_ident@sample)
   
   scTypeEval@dissimilarity[[method]] <- rr
      
   
   return(scTypeEval)
   
}



#' @title Run Consistency Analysis on Single-Cell Data
#' @description Computes internal validation consistency metrics for cell type annotations
#' in single-cell RNA-seq datasets. This function evaluates the robustness of annotations
#' using various distance metrics, normalization methods, and gene sets.
#'
#' @param scTypeEval An scTypeEval object containing single-cell expression data, metadata, and gene lists.
#' @param ident Character. Name of the column in `scTypeEval@metadata` containing the cell type labels.
#'   If `NULL`, defaults to `scTypeEval@active.ident`.
#' @param sample Character. Name of the column in `scTypeEval@metadata` containing sample identifiers.
#'   Required for pseudobulk analyses.
#' @param normalization.method Character. Method for normalizing the expression data. See \link[scTypeEval]{add.HVG} for more details.
#'   Options: `"Log1p"`, `"CLR"`, `"pearson"`. Default: `"Log1p"`. 
#' @param gene.list List. A list of gene sets to compute consistency metrics on.
#'    By default, each of the gene lists stored in `scTypeEval` are used recursively to compute consitency metrics.
#' @param pca Logical. Whether to perform PCA before computing metrics. Default: `FALSE`. This parameter will be turned to `TRUE` for GloScope data.type.
#'    `FALSE` will build consistency metrics directly on the genes within `gene.list`, while `TRUE` will do it on their principal components.
#' @param ndim Integer. Number of PCA dimensions to use to compute metrics if `pca = TRUE`. Default: `30`.
#' @param distance.method Character. Distance metric to use. Must be one of the predefined methods. Default: `"euclidean"`.
#'   Supported options:
#'   - `"euclidean"`: Euclidean distance, commonly used for measuring dissimilarity in high-dimensional spaces.
#'   - `"EMD"`: Earth Mover's Distance, computationally expensive but useful for comparing distributions (only recommended for pseudobulk data).
#'   - `"maximum"`: Chebyshev distance, considers the largest absolute difference across dimensions.
#'   - `"manhattan"`: Sum of absolute differences, often preferred when dealing with sparse data.
#'   - `"canberra"`: Weighted Manhattan distance, emphasizing smaller values.
#'   - `"binary"`: Simple matching coefficient for binary presence/absence data.
#'   - `"minkowski"`: Generalized distance metric (Manhattan and Euclidean are special cases).
#'   - `"Jaccard"`: Measures dissimilarity between binary vectors.
#'   - `"Weighted_Jaccard"`: Variant of Jaccard similarity that accounts for weighting.
#'   - `"gower"`: Distance measure handling mixed binary raw + normalized matrices.
#'   - `"bray-curtis"`: Measures dissimilarity between samples based on abundance.
#'   - `"cosine"`: Measures the cosine of the angle between two vectors.
#'   - `"pearson"`: Pearson correlation distance (1 - correlation coefficient).
#' @param IntVal.metric Character vector. Internal validation metrics to compute. By default, the following metrics will be run:
#'   `"silhouette"`, `"NeighborhoodPurity"`, `"ward.PropMatch"`, `"Leiden.PropMatch"`, `"modularity"`.
#'     
#'   All available options:
#'   - `"silhouette"`: Contrast intra-cluster tightness with inter-cluster separation. (Distance-based, Mean silhouette width per cell type)
#'   - `"modularity"`: Strength of intra-community density compared to random expectation in a network. (Graph structure, Global modularity score and per-cell-type modularity contribution)
#'   - `"ward.PropMatch"`: Normalized proportion of reference labels in the dominant cluster by Ward (Clustering-based, Proportion match of dominant labels within hierarchical clusters)
#'   - `"ward.NMI"`: Normalized shared information between Ward clusters and reference labels. (Global level metric)
#'   - `"ward.ARI"`: Adjusted Rand Index for Ward clustering. (Global level metric)
#'   - `"Leiden.PropMatch"`: Normalized proportion of reference labels in the dominant Leiden clustering. (Clustering-based, Local agreement)
#'   - `"Leiden.NMI"`: Normalized Mutual Information for Leiden partitions. (Global level metric)
#'   - `"Leiden.ARI"`: Adjusted Rand Index for Leiden clustering. (Global level metric)
#'   - `"NeighborhoodPurity"`: Proportion of K-Nearest Neighbors sharing the same label. (Local agreement, Mean proportion of neighbors sharing same identity)
#'   - `"GraphConnectivity"`: Proportion of elements within a group that belong to its largest connected component. (Graph structure, Size of the largest connected component relative to reference group size)
#'   - `"Orbital.centroid"`: Proportion of cells closer to their own cluster centroid than to any cluster other centroid. (Distance-based euclidean, Local agreement using Euclidean distance to centroids)
#'   - `"Orbita.medoid"`: Proportion of cells closer to their own cluster medoid than to any other cluster medoid. (Distance-based euclidean, Local agreement using Euclidean distances to medoid cell/sample)
#'   - `"inertia"`: Sum of squared distances from points to their assigned cluster centroid. (Distance-based)  
#'   - `"Xie-Beni"`: Ratio of intra-cluster dispersion to inter-cluster separation. (Distance-based)  
#'   - `"S.Dbw"`: Internal cluster validity index combining density and separation. (Distance-based)  
#'   - `"I"`: Trade-off between separation and cohesion in clusters. (Distance-based) 
#' @param data.type Character. Type of data input to perform.
#'   Options, one of:
#'   - `"sc"`: Single-cell data, where each cell is treated as an individual observation.
#'   - `"pseudobulk"`: Aggregated expression values per cell type and sample, useful to capture inter-sample variability.
#'   - `"pseudobulk_1vsall"`: Pseudobulk comparisons where each cell type is analyzed against all others, useful for detecting cell-type-specific consistency.
#'   - `"GloScope"`: Estimates the gene expression density from a low dimensional embedding of the UMI count measurements and measure its Kullback_Leibler (KL) divergency using \link{GloScope}.
#'   Default: `"sc"`.
#' @param min.samples Integer. Minimum number of samples required for intersample comparisons. Default: `5`.
#' @param min.cells Integer. Minimum number of cells required per population. Default: `10`.
#' @param KNNGraph_k Integer. Number of neighbors to consider in k-NN graph calculations. Default: `5`.
#' @param black.list Character vector. Genes to exclude from consistency analysis. Defaults to `scTypeEval@black.list`.
#' @param ncores Integer. Number of CPU cores to use for parallel processing. Default: `1`.
#' @param bparam Parallel backend parameter object for BiocParallel. If provided, overrides `ncores`.
#' @param progressbar Logical. Whether to display a progress bar during computation. Default: `TRUE`.
#' @param verbose Logical. Whether to print messages during execution. Default: `TRUE`.
#'
#' @return An updated `scTypeEval` object with consistency metrics stored in \code{consistency} slot.
#'
#' 
#' @details The function evaluates the consistency of cell type annotations using internal validation metrics.  
#' If `sample` is provided, it performs inter-sample consistency analysis (pseudobulk modes require multiple samples).  
#' The function supports different validation metrics, including silhouette score,  
#' network modularity, and nearest-neighbor-based purity measures. Users can specify a gene set  
#' for evaluation and optionally apply PCA before computing consistency metrics.  
#' Parallel processing is enabled when `ncores > 1`, using BiocParallel for efficiency. 
#' 
#' 
#' @examples
#' \dontrun{
#' # Run consistency analysis on single-cell data
#' sceval <- Run.Consistency(
#'   scTypeEval = sceval,
#'   ident = "cell_type",
#'   sample = "sample_id",
#'   normalization.method = "Log1p",
#'   IntVal.metric = c("silhouette", "modularity"),
#'   data.type = "pseudobulk",
#'   ncores = 2
#' )
#' }
#' 
#' @importClassesFrom Matrix dgCMatrix
#'
#' @export get.Consistency


get.Consistency <- function(scTypeEval,
                            dissimilarity.slot, 
                            IntVal.metric = c("silhouette",
                                              "NeighborhoodPurity",
                                              "ward.PropMatch",
                                              "Orbital.medoid",
                                              "Average.similarity"),
                            KNNGraph_k = 5,
                            verbose = TRUE
                            
){
   
   diss.assays <- .check_dissimilarityAssays(scTypeEval, slot = dissimilarity.slot)
   
   consist.list <- lapply(diss.assays,
                          function(da){
                             assay <- scTypeEval@dissimilarity[[da]]
                             
                             dist <- assay@dissimilarity
                             ident <- unlist(assay@ident)
                             
                             # compute internal validation metrics
                             if(verbose){message("Computing internal validation metrics ", da, " ... \n")}
                             con <- calculate_IntVal_metric(dist = dist,
                                                            metrics = IntVal.metric,
                                                            ident = ident,
                                                            KNNGraph_k = KNNGraph_k,
                                                            verbose = verbose)
                             dfl <- lapply(names(con),
                                          function(int){
                                             # build dissimilarity object
                                             r <- data.frame(measure = con[[int]],
                                                             consistency.metric = int)
                                             return(r)
                                          })
                             
                             df <- do.call(rbind, dfl) |>
                                dplyr::mutate(dissimilarity_method = assay@method,
                                              ident = names(assay@ident))
                             
                             return(df)
                          })
   
   consist <- do.call(rbind, consist.list) 
   
   return(consist)
   
}





#' Plot PCA Results from scTypeEval Object
#'
#' This function visualizes Principal Component Analysis (PCA) results stored in the \code{reductions} slot of an scTypeEval object.
#'
#' @param scTypeEval An \code{scTypeEval} object containing PCA results in the \code{reductions} slot.
#' @param gene.list Character vector. gene.list names to filter PCA plots. If \code{NULL}, all available PCA results are plotted.
#' @param data.type Character vector. Specifies which data types to include in the plots (e.g., \code{"sc"}, \code{"pseudobulk"}).
#' If \code{NULL}, all data types are included.
#' @param dims Integer vector of length 2. The principal component (PC) dimensions to plot (default: \code{c(1,2)}).
#' @param label Logical. Whether to add labels to the PCA plot (default: \code{TRUE}).
#' @param show.legend Logical. Whether to display a legend (default: \code{FALSE}).
#'
#' @return A named list of PCA plots (\link{ggplot2} objects) corresponding to different PCA analyses stored in
#'  \code{reductions} slot of an scTypeEval object by \link{add.PCA}.
#'
#' @seealso \link{add.PCA}
#' 
#' @examples
#' \dontrun{
#' # plot all PCA dimensional reductions
#' pca_plots <- plot.PCA(sceval) # scTypeEval object with previously run `add.PCA()`
#' 
#' # plot only some PCAs
#' pca_plots_filtered <- plot.PCA(sceval, # scTypeEval object with previously run `add.PCA()`
#'                               gene.list = "HVG", # plot only PCA for highly variable genes
#'                               data.type = "pseudobulk" # plot only PCA for pseudobulk data
#'                               ) 
#' }
#'
#' @export plot.PCA


plot.PCA <- function(scTypeEval,
                     gene.list = NULL,
                     data.type = NULL,
                     dims = c(1,2),
                     label = TRUE,
                     show.legend = FALSE) {
   
   # Ensure the reductions slot is not empty
   if (length(scTypeEval@reductions) == 0) {
      stop("The 'reductions' slot in the 'scTypeEval' object is empty.")
   }
   
   assays <- unlist(scTypeEval@reductions)
   pls <- list()  # Initialize an empty list to store plots
   
   # Loop through each ConsistencyAssay object in scTypeEval@reductions
   for (a in names(assays)) {
      assay <- assays[[a]]
      
      # Check for class validity
      if (!inherits(assay, "DimRed")) {
         stop(paste("Invalid object in reductions slot"))
      }
      
      # Apply filtering
      if (!is.null(gene.list) && length(intersect(gene.list, assay@gene.list)) == 0) next
      if (!is.null(data.type) && !assay@data.type %in% data.type) next
      if (assay@key != "PCA") next
      
      df <- assay@embeddings[, dims] |>
         as.data.frame() |>
         dplyr::mutate(ident = assay@ident)
      
      # compute variance of PCs
      vrs <- var_PCA(assay@embeddings)[dims]
      vrs <- round(vrs*100, 2)
      
      labs <- paste0("PC", dims, " (", vrs, "%)")
      
      pl <- helper.plot.PCA(df,
                            show.legend = show.legend,
                            label = label) +
         ggplot2::labs(x = labs[1],
                       y = labs[2],
                       title = a)
      
      pls[[a]] <- pl  # Store plot in the list
   }
   
   if(length(pls) == 0){
      warning("No PCA reductions found, run add.PCA before.\n")
   }
   
   return(pls)
}

#' Perform Hierarchical Clustering on scRNA-seq Data
#'
#' This function performs hierarchical clustering on single-cell RNA sequencing (scRNA-seq) data 
#' using specified normalization, distance, and clustering methods. It allows for analysis at 
#' both the single-cell and pseudobulk levels.
#'
#' @param scTypeEval An object containing the scRNA-seq data, metadata, and gene lists.
#' @param ident A character string specifying the metadata column used to group cells (e.g., cell type or annotation).
#' @param sample An optional character string specifying the metadata column for sample identification (required for pseudobulk analysis).
#' @param normalization.method A character vector specifying the normalization method. Options: `"Log1p"`, `"CLR"`, `"pearson"` (default: `"Log1p"`).
#' @param gene.list A named list of gene sets to be used for clustering. If `NULL`, gene lists stored in `scTypeEval` are used.
#' @param pca Logical; if `TRUE`, performs PCA dimensionality reduction before clustering (default: `FALSE`).
#' @param ndim Integer; number of principal components to retain if `pca = TRUE` (default: `30`).
#' @param distance.method A character string specifying the distance metric for clustering (default: `"euclidean"`).
#' @param hierarchy.method A character string specifying the hierarchical clustering method. Default is `"ward.D2"`. For more option see \link[stats]{hclust}.
#' @param data.type A character string indicating whether to use `"pseudobulk"` (aggregated samples), `"GloScope`, or `"sc"` (single-cell) data (default: `"pseudobulk"`).
#' @param min.samples Integer; minimum number of samples required for pseudobulk analysis (default: `5`).
#' @param min.cells Integer; minimum number of cells required per group for inclusion in the analysis (default: `10`).
#' @param black.list A vector of gene names to exclude from clustering. If `NULL`, uses the blacklist stored in `scTypeEval`.
#' @param ncores Integer; number of CPU cores to use for parallel processing (default: `1`).
#' @param bparam An optional `BiocParallelParam` object for controlling parallel execution. If provided, overrides ncores.
#' @param progressbar Logical; if `TRUE`, displays a progress bar during computation (default: `TRUE`).
#' @param verbose Logical; if `TRUE`, prints messages about progress and warnings (default: `TRUE`).
#'
#' @return A list of hierarchical clustering results, where each element corresponds to a different gene list.
#'
#' @details The function retrieves expression data for specified gene lists, applies normalization, 
#' and performs hierarchical clustering using the chosen distance metric and clustering method.
#' If `pca = TRUE`, PCA is performed before clustering to reduce dimensionality.
#'
#' @examples
#' \dontrun{
#' result <- get.hierarchy(scTypeEval = sceval, 
#'                         ident = "celltype",
#'                         sample = "sample_id",
#'                         normalization.method = "Log1p",
#'                         data.type = "pseudobulk",
#'                         distance.method = "euclidean",
#'                         hierarchy.method = "ward.D2")
#' # explore results
#' clusters <- stats::cutree(result[[HVG]], k = length(cell_types))
#' }
#'
#' @export get.hierarchy

get.hierarchy <- function(scTypeEval,
                          ident,
                          sample,
                          normalization.method = c("Log1p", "CLR", "pearson"),
                          gene.list = NULL,
                          pca = FALSE,
                          ndim = 30,
                          distance.method = "euclidean",
                          hierarchy.method = "ward.D2",
                          data.type = c("pseudobulk", "sc", "GloScope"),
                          min.samples = 5,
                          min.cells = 10,
                          black.list = NULL,
                          ncores = 1,
                          bparam = NULL,
                          progressbar = FALSE,
                          verbose = TRUE
){
   if(is.null(ident)){
      ident <- scTypeEval@active.ident
   }
   
   if(!ident %in% names(scTypeEval@metadata)){
      stop("Please provide a ident, i.e. a cell type or annotation to group cells included in metadata")
   }
   
   data.type = data.type[1]
   
   # retrieve ident and convert to factor
   ident.name <- ident
   ident <- scTypeEval@metadata[[ident]]
   ident <- purge_label(ident)
   ident <- factor(ident)
   
   if(length(unique(ident)) < 2){
      stop("Less than 2 cell type detected, Consistency metrics not possible")
   }
   
   if(data.type == "sc"){
      sample <- NULL
   }
   
   if(!is.null(sample)){
      if(!sample %in% names(scTypeEval@metadata)){
         stop("`sample` parameter not found in metadata colnames.")
      }
      # retrieve sample and convert to factor
      sample.name <- sample
      sample <- scTypeEval@metadata[[sample]]
      sample <- purge_label(sample)
      sample <- factor(sample)
      
      if(min.samples < 2){
         stop("For intersample comparison a minimum threshold of cell population
              in at least 2 samples is required, but more is recommended.")
      }
      
   } else {
      
      if(data.type != "sc"){
         stop("For pseudobulk or GloScope provide a dataset with multiple samples,
              and specifiy their respective column metadata in `sample` parameter")
      } else {
         if(verbose){message("Using dataset as a unique sample, computing consistency across cells.\n")}
      }
      sample.name <- NULL
   }
   
   
   if(!distance.method %in% distance_methods){
      stop(distance.method, " distance method not supported. Pick up one of: ", 
           paste(distance_methods, collapse = ", "))
   }
   # only run EMD on pseudobulk data.type
   if(distance.method == "EMD" && data.type == "sc"){
      warning("Earth mover distance (EMD) for single-cell (sc) data.type is
              highly computially expensive... not recommended, it will take long")
   }
   
   if(!data.type %in% data_type){
      stop(data.type, " data type conversion method not supported. Pick up one of: ", 
           paste(data_type, collapse = ", "))
   }
   
   if(data.type == "GloScope" && !pca){
      warning("data.type GloScope only run in PCA space, switching to pca=TRUE\n")
      pca <- TRUE
   }
   
   if(data.type == "GloScope" && !distance.method %in% Gloscope_dists){
      stop("data.type GloScope only possible with ",
           paste(Gloscope_dists, collapse = ", "),
           " distances methods.")
   }
   
   if(data.type != "GloScope" && distance.method %in% Gloscope_dists){
      stop(paste(Gloscope_dists, collapse = ", "),
           " only supported with GloScope data.type.")
   }
   
   # set normalization method
   normalization.method <- normalization.method[1]
   
   # set gene lists
   if(is.null(gene.list)){
      gene.list <- scTypeEval@gene.lists
      if(length(gene.list) == 0){
         stop("No gene list found to run consistency metrics.\n
              Add custom gene list or compute highly variable genes with add.HVG()\n")
      }
   } else {
      if(!all(names(gene.list) %in% names(scTypeEval@gene.lists))){
         stop("Some gene list names not included in scTypeEval object")
      }
      gene.list <- scTypeEval@gene.lists[gene.list]
   }
   
   if(is.null(black.list)){
      black.list <- scTypeEval@black.list
   }
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   # loop over each gene.list
   hier.list <- BiocParallel::bplapply(names(gene.list),
                                       BPPARAM = param,
                                       function(t){
                                          
                                          # get matrix
                                          keep <- rownames(scTypeEval@counts) %in% gene.list[[t]]
                                          mat <- scTypeEval@counts[keep,]
                                          
                                          # remove black list genes
                                          mat <- mat[!rownames(mat) %in% black.list,]
                                          
                                          hier <- hierarchy.helper(mat,
                                                                   ident = ident,
                                                                   sample = sample,
                                                                   normalization.method = normalization.method,
                                                                   distance.method = distance.method,
                                                                   hierarchy.method = hierarchy.method,
                                                                   data.type = data.type,
                                                                   pca = pca,
                                                                   ndim = ndim,
                                                                   min.samples = min.samples,
                                                                   min.cells = min.cells,
                                                                   verbose = verbose)
                                          return(hier)
                                       })
   
   if(pca){
      names(hier.list) <- paste(names(gene.list), "PCA", sep = ".")
   } else {
      names(hier.list) <- names(gene.list)
   }
   
   return(hier.list)
}

#' Compute K-Nearest Neighbors (KNN) Graph on scRNA-seq Data
#'
#' This function computes a KNN graph based on single-cell RNA sequencing (scRNA-seq) data,
#' allowing for consistency analysis within and across cell types. The analysis can be performed
#' on single-cell or pseudobulk data.
#'
#' @param scTypeEval An object containing scRNA-seq expression data, metadata, and gene lists.
#' @param ident A character string specifying the metadata column used to group cells (e.g., cell type or annotation).
#' @param sample An optional character string specifying the metadata column for sample identification (required for pseudobulk analysis).
#' @param normalization.method A character vector specifying the normalization method. Options: `"Log1p"`, `"CLR"`, `"pearson"` (default: `"Log1p"`).
#' @param gene.list A named list of gene sets to be used for clustering. If `NULL`, the gene lists stored in `scTypeEval` are used.
#' @param pca Logical; if `TRUE`, performs PCA dimensionality reduction before computing distances (default: `FALSE`).
#' @param ndim Integer; number of principal components to retain if `pca = TRUE` (default: `30`).
#' @param distance.method A character string specifying the distance metric for computing the KNN graph (default: `"euclidean"`).
#' @param KNNGraph_k Integer; the number of neighbors to consider in the KNN graph (default: `5`).
#' @param data.type A character string indicating whether to use `"pseudobulk"` (aggregated samples), `"GloScope`, or `"sc"` (single-cell) data (default: `"pseudobulk"`).
#' @param min.samples Integer; minimum number of samples required for pseudobulk analysis (default: `5`).
#' @param min.cells Integer; minimum number of cells required per group for inclusion in the analysis (default: `10`).
#' @param black.list A vector of gene names to exclude from analysis. If `NULL`, uses the blacklist stored in `scTypeEval`.
#' @param ncores Integer; number of CPU cores to use for parallel processing (default: `1`).
#' @param bparam An optional `BiocParallelParam` object for controlling parallel execution.
#' @param progressbar Logical; if `TRUE`, displays a progress bar during computation (default: `TRUE`).
#' @param verbose Logical; if `TRUE`, prints messages about progress and warnings (default: `TRUE`).
#'
#' @return A list containing KNN graph-based consistency metrics for each gene list.
#'
#' @details This function constructs a KNN graph from scRNA-seq expression data using a chosen
#' distance metric. It normalizes gene expression, applies optional PCA, and computes KNN relationships
#' based on cell or pseudobulk similarities. The function assesses cell type consistency within the graph.
#'
#' @examples
#' 
#' \dontrun{
#' result <- get.NN(scTypeEval = sceval, 
#'                  ident = "celltype",
#'                  sample = "sample_id",
#'                  normalization.method = "Log1p",
#'                  data.type = "pseudobulk",
#'                  distance.method = "euclidean",
#'                  KNNGraph_k = 5)
#' }
#'
#' @export get.NN

get.NN <- function(scTypeEval,
                   ident,
                   sample,
                   normalization.method = c("Log1p", "CLR", "pearson"),
                   gene.list = NULL,
                   pca = FALSE,
                   ndim = 30,
                   distance.method = "euclidean",
                   KNNGraph_k = 5,
                   data.type = c("pseudobulk", "sc", "GloScope"),
                   min.samples = 5,
                   min.cells = 10,
                   black.list = NULL,
                   ncores = 1,
                   bparam = NULL,
                   progressbar = FALSE,
                   verbose = TRUE
){
   if(is.null(ident)){
      ident <- scTypeEval@active.ident
   }
   
   if(!ident %in% names(scTypeEval@metadata)){
      stop("Please provide a ident, i.e. a cell type or annotation to group cells included in metadata")
   }
   
   data.type = data.type[1]
   
   # retrieve ident and convert to factor
   ident.name <- ident
   ident <- scTypeEval@metadata[[ident]]
   ident <- purge_label(ident)
   ident <- factor(ident)
   
   if(length(unique(ident)) < 2){
      stop("Less than 2 cell type detected, Consistency metrics not possible")
   }
   
   if(data.type == "sc"){
      sample <- NULL
   }
   
   if(!is.null(sample)){
      if(!sample %in% names(scTypeEval@metadata)){
         stop("`sample` parameter not found in metadata colnames.")
      }
      # retrieve sample and convert to factor
      sample.name <- sample
      sample <- scTypeEval@metadata[[sample]]
      sample <- purge_label(sample)
      sample <- factor(sample)
      
      if(min.samples < 2){
         stop("For intersample comparison a minimum threshold of cell population
              in at least 2 samples is required, but more is recommended.")
      }
      
   } else {
      
      if(data.type != "sc"){
         stop("For pseudobulk or GloScope provide a dataset with multiple samples,
              and specifiy their respective column metadata in `sample` parameter")
      } else {
         if(verbose){message("Using dataset as a unique sample, computing consistency across cells.\n")}
      }
      sample.name <- NULL
   }
   
   
   if(!distance.method %in% distance_methods){
      stop(distance.method, " distance method not supported. Pick up one of: ", 
           paste(distance_methods, collapse = ", "))
   }
   # only run EMD on pseudobulk data.type
   if(distance.method == "EMD" && data.type == "sc"){
      warning("Earth mover distance (EMD) for single-cell (sc) data.type is
              highly computially expensive... not recommended, it will take long")
   }
   
   if(!data.type %in% data_type){
      stop(data.type, " data type conversion method not supported. Pick up one of: ", 
           paste(data_type, collapse = ", "))
   }
   
   if(data.type == "GloScope" && !pca){
      warning("data.type GloScope only run in PCA space, switching to pca=TRUE\n")
      pca <- TRUE
   }
   
   if(data.type == "GloScope" && !distance.method %in% Gloscope_dists){
      stop("data.type GloScope only possible with ",
           paste(Gloscope_dists, collapse = ", "),
           " distances methods.")
   }
   
   if(data.type != "GloScope" && distance.method %in% Gloscope_dists){
      stop(paste(Gloscope_dists, collapse = ", "),
           " only supported with GloScope data.type.")
   }
   
   # set normalization method
   normalization.method <- normalization.method[1]
   
   # set gene lists
   if(is.null(gene.list)){
      gene.list <- scTypeEval@gene.lists
      if(length(gene.list) == 0){
         stop("No gene list found to run consistency metrics.\n
              Add custom gene list or compute highly variable genes with add.HVG()\n")
      }
   } else {
      if(!all(names(gene.list) %in% names(scTypeEval@gene.lists))){
         stop("Some gene list names not included in scTypeEval object")
      }
      gene.list <- scTypeEval@gene.lists[gene.list]
   }
   
   if(is.null(black.list)){
      black.list <- scTypeEval@black.list
   }
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   # loop over each gene.list
   nn.list <- BiocParallel::bplapply(names(gene.list),
                                     BPPARAM = param,
                                     function(t){
                                        
                                        # get matrix
                                        keep <- rownames(scTypeEval@counts) %in% gene.list[[t]]
                                        mat <- scTypeEval@counts[keep,]
                                        
                                        # remove black list genes
                                        mat <- mat[!rownames(mat) %in% black.list,]
                                        
                                        nn <- nn.helper(mat,
                                                        ident = ident,
                                                        sample = sample,
                                                        normalization.method = normalization.method,
                                                        distance.method = distance.method,
                                                        KNNGraph_k = KNNGraph_k,
                                                        data.type = data.type,
                                                        pca = pca,
                                                        ndim = ndim,
                                                        min.samples = min.samples,
                                                        min.cells = min.cells,
                                                        verbose = verbose)
                                        return(nn)
                                     })
   
   if(pca){
      names(nn.list) <- paste(names(gene.list), "PCA", sep = ".")
   } else {
      names(nn.list) <- names(gene.list)
   }
   
   return(nn.list)
}



dx.BestHit <- function(scTypeEval,
                       ident = NULL,
                       method = c("Mutual.Score", "Mutual.Match"),
                       sample = NULL,
                       data.type = c("pseudobulk, sc"),
                       gene.list = NULL,
                       black.list = NULL,
                       min.cells = 10,
                       min.samples = 5,
                       ncores = 1,
                       bparam = NULL,
                       progressbar = FALSE,
                       verbose = TRUE
                       
){
   
   data.type <- data.type[1]
   if(!data.type %in% data_type){
      stop(data.type, " data type conversion method not supported. Pick up one of: ", 
           paste(data_type[1:2], collapse = ", "))
   } else if(data.type %in% data_type[3:4]){
      stop(paste(data_type[3:4], collapse = ", "), "not supported for mutual BestHit.") 
   }
   
   if(is.null(ident)){
      ident <- scTypeEval@active.ident
   }
   
   if(!ident %in% names(scTypeEval@metadata)){
      stop("Please provide a ident, i.e. a cell type or annotation to group cells included in metadata")
   }
   
   # retrieve ident and convert to factor
   ident.name <- ident
   ident <- scTypeEval@metadata[[ident]]
   ident <- purge_label(ident)
   ident <- factor(ident)
   
   if(length(unique(ident)) < 2){
      stop("Less than 2 cell type detected, Consistency metrics not possible")
   }
   
   if(!is.null(sample)){
      if(!sample %in% names(scTypeEval@metadata)){
         stop("`sample` parameter not found in metadata colnames.")
      }
      # retrieve sample and convert to factor
      sample.name <- sample
      sample <- scTypeEval@metadata[[sample]]
      sample <- purge_label(sample)
      sample <- factor(sample)
   } else {
      stop("BestHit consistency requires multiple samples and specify it in `sample` parameter.")
   }
   
   # set gene lists
   if(is.null(gene.list)){
      gene.list <- scTypeEval@gene.lists
      if(length(gene.list) == 0){
         stop("No gene list found to run consistency metrics.\n
              Add custom gene list or compute highly variable genes with add.HVG()\n")
      }
   } else {
      if(!all(names(gene.list) %in% names(scTypeEval@gene.lists))){
         stop("Some gene list names not included in scTypeEval object")
      }
      gene.list <- scTypeEval@gene.lists[gene.list]
   }
   
   if(is.null(black.list)){
      black.list <- scTypeEval@black.list
   }
   
   # set methods
   if(any(!method %in% mutual_method)){
      stop("Not supported consistency metrics, please provide both or either: ", mutual_method)
   }
   
   if(data.type == "sc" && "Mutual.Match" %in% method){
      warning("Mutual.Match consistency only supported for pseudobulk data.type, not running")
      method <- method[method != "Mutual.Match"]
   }
   
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   # loop over each gene.list
   dx.list <- lapply(names(gene.list),
                     function(t){
                        
                        # get matrix
                        keep <- rownames(scTypeEval@counts) %in% gene.list[[t]]
                        mat <- scTypeEval@counts[keep,]
                        
                        # remove black list genes
                        mat <- mat[!rownames(mat) %in% black.list,]
                        
                        con <- dx.bestHit(mat = mat,
                                          ident = ident,
                                          sample = sample,
                                          data.type = data.type,
                                          method = method,
                                          min.cells = min.cells,
                                          min.samples = min.samples,
                                          bparam = param)
                        return(con)
                     })
   
   names(dx.list) <- names(gene.list)
   
   
   return(dx.list)
}




#' Plot PCA Results from scTypeEval Object
#'
#' This function visualizes Principal Component Analysis (PCA) results stored in the \code{reductions} slot of an scTypeEval object.
#'
#' @param scTypeEval An \code{scTypeEval} object containing PCA results in the \code{reductions} slot.
#' @param gene.list Character vector. gene.list names to filter PCA plots. If \code{NULL}, all available PCA results are plotted.
#' @param data.type Character vector. Specifies which data types to include in the plots (e.g., \code{"sc"}, \code{"pseudobulk"}).
#' If \code{NULL}, all data types are included.
#' @param dims Integer vector of length 2. The principal component (PC) dimensions to plot (default: \code{c(1,2)}).
#' @param label Logical. Whether to add labels to the PCA plot (default: \code{TRUE}).
#' @param show.legend Logical. Whether to display a legend (default: \code{FALSE}).
#'
#' @return A named list of PCA plots (\link{ggplot2} objects) corresponding to different PCA analyses stored in
#'  \code{reductions} slot of an scTypeEval object by \link{add.PCA}.
#'
#' @seealso \link{add.MDS}
#' 
#' @examples
#' \dontrun{
#' # plot all PCA dimensional reductions
#' pca_plots <- plot.MDS(sceval) # scTypeEval object with previously run `add.PCA()`
#' 
#' # plot only some PCAs
#' mds_plots_filtered <- plot.MDS(sceval, # scTypeEval object with previously run `add.PCA()`
#'                               gene.list = "HVG", # plot only PCA for highly variable genes
#'                               data.type = "pseudobulk" # plot only PCA for pseudobulk data
#'                               ) 
#' }
#'
#' @export plot.MDS


plot.MDS <- function(scTypeEval,
                     dissimilarity.slot = "all",
                     label = TRUE,
                     dims = c(1,2),
                     show.legend = FALSE) {
   
   diss.assays <- .check_dissimilarityAssays(scTypeEval, slot = dissimilarity.slot)
   
   pls <- lapply(diss.assays,
                 function(da){
                    assay <- scTypeEval@dissimilarity[[da]]
                    
                    dist <- assay@dissimilarity
                    ident <- unlist(assay@ident)
                    title <- paste(da, names(assay@ident), sep = " - ")
                    
                    mds <- stats::cmdscale(dist,
                                           k = max(dims),
                                           list. = FALSE)
                    
                    df <- mds[, dims] |>
                       as.data.frame() |>
                       dplyr::mutate(ident = ident)
                    
                    labs <- paste0("Dim", dims)
                    
                    pl <- helper.plot.scatter(df,
                                              show.legend = show.legend,
                                              label = label) +
                       ggplot2::labs(x = labs[1],
                                     y = labs[2],
                                     title = title)
                    return(pl)
                 })
   
   names(pls) <- diss.assays
   
   if(length(pls) == 1){
      pls <- unlist(pls)
      }
   
   return(pls)
}


plot.Heatmap <- function(scTypeEval,
                         dissimilarity.slot = "all",
                         ... # other params for pheatmap
) {
   
   if (!requireNamespace("pheatmap", quietly = TRUE)) {
      message("Installing required package pheatmap...\n")
      install.packages("pheatmap")
   }
   
   diss.assays <- .check_dissimilarityAssays(scTypeEval, slot = dissimilarity.slot)
   
   pls <- lapply(diss.assays,
                 function(da){
                    assay <- scTypeEval@dissimilarity[[da]]
                    title <- paste(da, names(assay@ident), sep = " - ")
                    
                    d_mat <- as.matrix(assay@dissimilarity)
                    
                    annot <- data.frame(id = colnames(d_mat)) %>% 
                       mutate(ident = factor(unlist(assay@ident)),
                              sample = assay@sample) %>% 
                       tibble::column_to_rownames("id") %>% 
                       arrange(ident) %>% 
                       select(ident)
                    
                    pmat <- d_mat[rownames(annot), rownames(annot)]
                    
                    gaps <- table(annot$ident) %>% as.vector()
                    gaps <- cumsum(gaps)
                  
                    # Create the heatmap
                    pheatmap::pheatmap(pmat, 
                                       scale = scale,
                                       cluster_rows = F,
                                       cluster_cols = F,
                                       annotation_row = annot,
                                       annotation_col = annot,
                                       show_rownames = F,
                                       show_colnames = F,
                                       gaps_row = gaps,
                                       gaps_col = gaps,
                                       main = title,
                                       ...
                    )
                    
                    return(pl)
                 })
   
   names(pls) <- diss.assays
   
   if(length(pls) == 1){
      pls <- unlist(pls)
   }
   
   return(pls)
}
