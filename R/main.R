# Create scTypeEval object

#' @title Create an scTypeEval object for evaluating cell type classification.
#'
#' @description This function initializes an `scTypeEval` object from various input formats, including Seurat, 
#' SingleCellExperiment, or raw count matrices. It ensures compatibility by validating input types 
#' and structures before constructing the object.
#'
#' @param matrix A Seurat object, SingleCellExperiment object, dense/sparse count
#'   matrix, or `NULL`. If `NULL`, an empty object is initialized.
#' @param metadata A metadata dataframe. Required if `matrix` is a raw matrix or
#'   `NULL`. Must have as many rows as the number of cells (columns of the count matrix).
#' @param gene.lists A named list of gene sets to use in the evaluation (default: empty list).
#' @param black.list A character vector of genes to exclude from analysis (default: empty).
#' @param active.ident The active identity class or cluster label (optional). If provided,
#'   it is validated against the metadata.
#'
#' @return An `scTypeEval` object containing:
#' \itemize{
#'   \item \code{counts}: A sparse count matrix (dgCMatrix).
#'   \item \code{metadata}: A dataframe with metadata for each cell.
#'   \item \code{gene.lists}: A list of gene sets used in classification.
#'   \item \code{black.list}: A vector of excluded genes.
#'   \item \code{active.ident}: Active cluster identity (if provided).
#' }
#'
#' @details
#' - If `matrix = NULL`, the function initializes an empty `dgCMatrix` with 0 rows and 0 columns,
#'   and requires a metadata dataframe.
#' - For Seurat and SingleCellExperiment objects, counts and metadata are extracted automatically.
#' - For raw matrices, metadata must be provided explicitly.
#' - The function validates that the number of metadata rows matches the number of cells (matrix columns).
#' - If `active.ident` is provided, it is checked for consistency with the metadata.
#'
#' @examples
#' \dontrun{
#' # From count matrix and metadata dataframe
#' sceval <- create.scTypeEval(matrix = count_matrix, metadata = metadata)
#'
#' # From Seurat object
#' sceval <- create.scTypeEval(seurat_obj)
#'
#' # From SingleCellExperiment object
#' sceval <- create.scTypeEval(sce_obj)
#'
#' # Initialize an empty object (requires metadata)
#' empty_eval <- create.scTypeEval(matrix = NULL, metadata = data.frame(celltype = character()))
#' }
#' 
#' @importClassesFrom Matrix dgCMatrix
#'
#' @export create.scTypeEval


create.scTypeEval <- function(matrix = NULL,
                              metadata = NULL,
                              gene.lists = list(),
                              black.list = character(),
                              active.ident = NULL) {
   
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
                                  version = as.character(packageVersion("scTypeEval"))
   )
   if(!is.null(active.ident)){
      tmp <- .check_ident(scTypeEval_obj,
                          active.ident,
                          verbose = FALSE)
   }
   
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


#' @title Process and normalize data within an scTypeEval object
#'
#' @description Runs processing on the count matrix stored in an `scTypeEval`
#' object by aggregating, filtering, and normalizing data. 
#' Results are stored as `DataAssay` objects within the `scTypeEval`.
#'
#' @param scTypeEval An `scTypeEval` object generated by \code{create.scTypeEval}.
#' @param ident A column name from metadata to use as the identity class (e.g. cell type, cluster). 
#'   If `NULL`, defaults to \code{scTypeEval@active.ident}.
#' @param sample A column name from metadata indicating sample identifiers, required.
#' @param normalization.method A string specifying the normalization method to apply 
#'   (default: `"Log1p"`).
#' @param aggregation Method to group cells, either `"single-cell"` or `"pseudobulk"`. Default is both.
#' @param min.samples Minimum number of samples required to retain a cell type (default: 5).
#' @param min.cells Minimum number of cells required to retain a cell type in a sample (default: 10).
#' @param verbose Logical indicating whether to print progress messages (default: TRUE).
#'
#' @return An updated `scTypeEval` object containing:
#' \itemize{
#'   \item \code{data}: A list of `DataAssay` objects, one for single-cell other for pseudobulk.
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates and sets the identity (\code{ident}) and sample grouping (\code{sample}).
#'   \item Iterates over predefined \code{aggregation_types}: single-cell and pseudobulk.
#'   \item Extracts and filters the count matrix using 
#'         \code{min.samples} and \code{min.cells} thresholds.
#'   \item Normalizes the resulting matrix.
#'   \item Wraps the processed data into `DataAssay` objects.
#'   \item Stores the list of processed assays inside the `scTypeEval` object.
#' }
#'
#' @examples
#' \dontrun{
#' # Run data processing on an scTypeEval object
#' sceval <- create.scTypeEval(seurat_obj)
#' sceval <- Run.ProcessingData(sceval, ident = "celltype", sample = "patient_id")
#' }
#'
#' @export Run.ProcessingData


Run.ProcessingData <- function(scTypeEval,
                               ident = NULL,
                               sample = NULL,
                               aggregation = c("single-cell", "pseudobulk"),
                               normalization.method = "Log1p",
                               min.samples = 5,
                               min.cells = 10,
                               verbose = TRUE){
   #### preflights checks
   ident.name <- ident
   if(is.null(ident.name)){
      ident.name <- scTypeEval@active.ident
   }
   ident <- .check_ident(scTypeEval, ident, verbose = verbose)
   sample <- .check_sample(scTypeEval, sample, verbose = verbose)
   
   # set normalization method
   normalization.method <- normalization.method[1]
   
   if(any(!tolower(aggregation) %in% aggregation_types)){
      stop("Invalid aggregation type. Valid options are: ",
           paste(aggregation_types, collapse = ", "))
   }
   
   #### RUN 
   data.list <- lapply(aggregation,
                       function(ag){
                          if(verbose){message("# Processing data for ", ag, " ... \n")}
                          # Transform data
                          if(verbose){message("   Transforming and filtering count matrix... \n")}
                          mat <- get.matrix(scTypeEval@counts,
                                            ident = ident,
                                            sample = sample,
                                            aggregation = ag,
                                            min.samples = min.samples,
                                            min.cells = min.cells)
                          
                          # Normalize data
                          if(verbose){message("   Normalizing count matrix via " , normalization.method, "... \n")}
                          norm.mat <- Normalize_data(mat@matrix,
                                                     method = normalization.method)
                          
                          rr <- methods::new("DataAssay",
                                             matrix = norm.mat,
                                             aggregation = ag,
                                             group = mat@group,
                                             ident = setNames(list(mat@ident), ident.name),
                                             sample = mat@sample)
                          return(rr)
                       })
   names(data.list) <- aggregation
   
   # add to scTypeEval object
   scTypeEval@data <- data.list
   
   return(scTypeEval)
   
}

#' @title Add externally processed data to an scTypeEval object
#'
#' @description Adds a user-supplied processed dataset (e.g. single-cell or pseudobulk)
#' into an `scTypeEval` object as a `DataAssay`. Supports filtering and consistency checks
#' on cell type and sample annotations.
#'
#' @param scTypeEval An `scTypeEval` object generated by \code{create.scTypeEval}.
#' @param data A count matrix (dense or sparse) containing processed expression values 
#'   (either single-cell or pseudobulk aggregated).
#' @param aggregation A string specifying the aggregation type. Must be one of the 
#'   supported: `"single-cell"` or `"pseudobulk"`.
#' @param ident A vector of cell identities corresponding to the columns of \code{data}. 
#'   Used to define grouping (e.g. cell type). Required.
#' @param ident.name A string specifying the name under which the provided \code{ident}
#'   will be stored (default: `"custom"`).
#' @param sample A vector of sample identifiers corresponding to the columns of \code{data}. 
#'   Required unless already encoded in the object.
#' @param filter Logical indicating whether to filter the data based on \code{min.samples}
#'   and \code{min.cells}. Only allowed when \code{aggregation = "single-cell"} (default: FALSE).
#' @param min.samples Minimum number of samples required to retain a feature (default: 5).
#' @param min.cells Minimum number of cells required to retain a feature (default: 10).
#' @param verbose Logical indicating whether to print progress messages (default: TRUE).
#'
#' @return An updated `scTypeEval` object containing:
#' \itemize{
#'   \item \code{data}: A new `DataAssay` object stored under the specified aggregation type.
#' }
#'
#' @details
#' - The function validates that identity and sample annotations match the dimensions of the input data.
#' - For \code{aggregation = "single-cell"}, optional filtering removes groups with too few 
#'   samples or cells.
#' - For \code{aggregation = "pseudobulk"}, the function checks that each (sample, identity) 
#'   pair occurs exactly once (i.e., fully aggregated).
#' - Processed data is wrapped in a `DataAssay` object and added to the `scTypeEval`.
#'
#' @examples
#' \dontrun{
#' # Add single-cell processed data
#' sceval <- create.scTypeEval(seurat_obj)
#' sceval <- Add.ProcessedData(
#'   sceval,
#'   data = expr_matrix,
#'   aggregation = "single-cell",
#'   ident = celltypes,
#'   sample = patients,
#'   filter = TRUE
#' )
#'
#' # Add pseudobulk data
#' sceval <- Add.ProcessedData(
#'   sceval,
#'   data = pseudobulk_matrix,
#'   aggregation = "pseudobulk",
#'   ident = bulk_idents,
#'   sample = bulk_samples
#' )
#' }
#' 
#' @seealso \link{Run.ProcessingData}
#'
#' @export Add.ProcessedData


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
                      matrix = data,
                      aggregation = aggregation,
                      group = group,
                      ident = setNames(list(ident), ident.name),
                      sample = sample)
   
   # add to scTypeEval object
   scTypeEval@data[[aggregation]]<- rr
   
   return(scTypeEval)
   
}

#' @title Identify and add highly variable genes (HVG) to an scTypeEval object
#'
#' @description Detects highly variable genes (HVGs) from the normalized
#' single-cell data stored in an `scTypeEval` object. The identified HVGs are
#' stored in the `gene.lists` slot under `"HVG"`.
#'
#' @param scTypeEval An `scTypeEval` object containing normalized data in the
#'   `"single-cell"` slot (see \code{Run.ProcessingData}).
#' @param var.method Character string specifying the method for identifying highly
#'   variable genes. Options: `"scran"` (default) or `"basic"`.
#' @param ngenes Integer specifying the number of highly variable genes to retain
#'   (default: `2000`).
#' @param sample Logical indicating whether to leverage sample information when
#'   computing HVGs. If `TRUE`, the \code{sample} annotation stored in data is used. If `FALSE`, HVGs are computed without sample grouping.
#' @param aggregation Method to group cells stored in `scTypeEval@data`, either `"single-cell"` or `"pseudobulk"`. Default is `"single-cell"`.
#' @param black.list A character vector of genes to exclude from HVG selection.
#'   If `NULL`, uses the object’s internal blacklist (`scTypeEval@black.list`).
#' @param ncores Integer specifying the number of CPU cores to use for parallel
#'   processing (default: `1`).
#' @param bparam A \code{BiocParallel} backend parameter object for parallelization.
#'   If provided, overrides \code{ncores}.
#' @param progressbar Logical, whether to display a progress bar during computation
#'   (default: `FALSE`).
#' @param verbose Logical, whether to print progress messages (default: `TRUE`).
#' @param ... Additional arguments passed to internal HVG computation functions.
#'
#' @return The modified `scTypeEval` object with HVGs added to
#'   \code{scTypeEval@gene.lists[["HVG"]]}.
#'
#' @details
#' - Requires that normalized single-cell data has been generated with
#'   \code{Run.ProcessingData}.
#' - Genes present in the blacklist (\code{black.list}) are removed before HVG selection.
#' - Available HVG methods:
#'   \itemize{
#'     \item \code{"scran"}: Uses the \pkg{scran} \link[scran]{modelGeneVar} function to
#'     model gene-specific variance and identify biologically variable genes.
#'     \item \code{"basic"}: A simple variance-to-mean approach ranking genes by
#'     coefficient of variation, selecting the top \code{ngenes}.
#'   }
#'   
#' @seealso \link{Run.ProcessingData}, \link{Add.ProcessedData}
#'
#' @examples
#' \dontrun{
#' # Create scTypeEval object and run preprocessing
#' sceval <- create.scTypeEval(seurat_obj)
#' sceval <- Run.ProcessingData(sceval)
#'
#' # Compute highly variable genes
#' sceval <- Run.HVG(sceval, var.method = "scran", ngenes = 2000)
#' }
#'
#' @export Run.HVG


Run.HVG <- function(scTypeEval,
                    var.method = "scran",
                    ngenes = 2000,
                    sample = TRUE,
                    aggregation = "single-cell",
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
   
   if(!tolower(aggregation) %in% aggregation_types){
      stop("Invalid aggregation type. Valid options are: ",
           paste(aggregation_types, collapse = ", "))
   }
   # normalized matrix
   norm.mat <- scTypeEval@data[[aggregation]]
   if(is.null(norm.mat)){
      stop("No normalization slot found. Please run before `Run.ProcessingData()`.\n")
   }
   mat <- norm.mat@matrix
   
   if(sample){
      sample <- norm.mat@sample
   } else {
      warning("Not leveraging sample information for computing HVG.\n")
      sample <- NULL
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
                                   bparam = param,
                                   ...),
                 "scran" = get.GeneVar(norm.mat = mat,
                                       sample = sample,
                                       ngenes = ngenes,
                                       bparam = param,
                                       ...),
                 stop(var.method, " not supported for getting variable genes.")
   )
   
   
   scTypeEval@gene.lists[["HVG"]] <- hgv
   
   return(scTypeEval)
}


#' @title Identify and add marker genes to an scTypeEval object
#'
#' @description Identifies cell type marker genes from normalized single-cell data
#' stored in an `scTypeEval` object. The identified markers are stored in the
#' `gene.lists` slot under the chosen method (e.g. `"scran.findMarkers"`).
#'
#' @param scTypeEval An `scTypeEval` object containing normalized data in the
#'   `"single-cell"` slot (see \code{Run.ProcessingData}).
#' @param method A character string specifying the marker gene identification
#'   method. Currently supported:
#'   \itemize{
#'     \item `"scran.findMarkers"` — Uses \pkg{scran}'s
#'     \link[scran]{findMarkers} for differential expression analysis.
#'   }
#'   Default: `"scran.findMarkers"`.
#' @param ngenes.celltype Integer specifying the max number of marker genes to retain
#'   per cell type (default: `50`).
#' @param aggregation Method to group cells stored in `scTypeEval@data`, either `"single-cell"` or `"pseudobulk"`. Default is `"single-cell"`.
#' @param black.list A character vector of genes to exclude from marker selection.
#'   If `NULL`, uses the object’s internal blacklist (`scTypeEval@black.list`).
#' @param ncores Integer specifying the number of cores to use for parallel
#'   processing (default: `1`).
#' @param bparam Optional. A \code{BiocParallel} parameter object for controlling
#'   parallel computation. If provided, overrides \code{ncores}.
#' @param progressbar Logical, whether to display a progress bar during computation
#'   (default: `FALSE`).
#' @param verbose Logical, whether to print messages during execution (default: `TRUE`).
#' @param ... Additional arguments passed to the underlying marker detection function.
#'
#' @return The modified `scTypeEval` object with marker genes added to
#'   \code{scTypeEval@gene.lists[[method]]}.
#'
#' @details
#' - Requires that normalized single-cell data has been generated with
#'   \code{Run.ProcessingData}.
#' - Both cell identities and sample annotations are automatically extracted
#'   from the normalized data.
#' - Genes present in the blacklist (\code{black.list}) are removed before marker selection.
#' - For `"scran.findMarkers"`, the \pkg{scran} method \code{findMarkers} is applied
#'   to identify differentially expressed genes per cell type while adjusting for sample effects.
#'   
#' @seealso \link{Run.ProcessingData}, \link{Add.ProcessedData}
#'
#' @examples
#' \dontrun{
#' # Create scTypeEval object and run preprocessing
#' sceval <- create.scTypeEval(seurat_obj)
#' sceval <- Run.ProcessingData(sceval)
#'
#' # Identify marker genes per cell type
#' sceval <- Run.GeneMarkers(sceval, method = "scran.findMarkers", ngenes.celltype = 50)
#' }
#'
#' @export Run.GeneMarkers



Run.GeneMarkers <- function(scTypeEval,
                            method = c("scran.findMarkers"),
                            ngenes.celltype = 50,
                            aggregation = "single-cell",
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
   
   if(!tolower(aggregation) %in% aggregation_types){
      stop("Invalid aggregation type. Valid options are: ",
           paste(aggregation_types, collapse = ", "))
   }
   # normalized matrix
   norm.mat <- scTypeEval@data[[aggregation]]
   if(is.null(norm.mat)){
      stop("No normalization slot found. Please run before `Run.ProcessingData()`.\n")
   }
   mat <- norm.mat@matrix
   ident <- norm.mat@ident[[1]]
   sample <- norm.mat@sample
   
   # remove blacked listed genes
   if(!is.null(black.list) && verbose){message("Filtering out black listed genes... \n")}
   mat <- mat[!rownames(mat) %in% black.list,]
   
   if(verbose){message("Computing cell type markers for ", names(norm.mat@ident),  "... \n")}
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





#' @title Perform PCA on Processed Data and Store Results in scTypeEval Object
#'
#' @description This function computes Principal Component Analysis (PCA) for each processed
#' data aggregation stored in the \code{scTypeEval} object (e.g., single-cell or pseudobulk).
#' The resulting PCA embeddings and loadings are stored in the \code{reductions} slot of the object.
#'
#' @param scTypeEval A \code{scTypeEval} object containing processed expression data. See \code{Run.ProcessingData} or \code{Add.ProcessedData}).
#' @param gene.list Named list of character vectors. Each element defines a gene set for PCA analysis.
#' If \code{NULL}, all pre-defined gene lists stored in \code{scTypeEval} are used.
#' @param black.list Character vector of genes to exclude from PCA. If \code{NULL},
#' the blacklist stored in \code{scTypeEval} is used.
#' @param ndim Integer. Number of principal components to compute (default: 30).
#' @param verbose Logical. Whether to print progress messages during computation (default: \code{TRUE}).
#'
#' @return The modified \code{scTypeEval} object with PCA results stored in the
#' \code{reductions} slot for each processed data aggregation.
#'
#' @details
#' This function runs PCA on all processed data slots within the \code{scTypeEval} object.
#' Each PCA result is stored as a \code{DimRed} assay containing:
#' \itemize{
#'   \item \code{embeddings}: PCA coordinates of samples/cells.
#'   \item \code{feature.loadings}: Loadings of features (genes) on each PC.
#'   \item \code{gene.list}: The gene sets used for PCA.
#'   \item \code{black.list}: The genes excluded from PCA.
#'   \item \code{aggregation}: The aggregation type (e.g., \code{"single-cell"}, \code{"pseudobulk"}).
#'   \item \code{ident}, \code{sample}, and \code{group}: Metadata carried over from processed data.
#' }
#'
#' @examples
#' \dontrun{
#' # Run PCA on all processed aggregations
#' sceval <- Run.PCA(sceval,
#'                   gene.list = list(mygenes = c("Cd4", "Cd8a", "Foxp3")),
#'                   ndim = 20)
#' }
#'
#' @seealso \link{Add.ProcessedData}
#'
#' @export Run.PCA



# obtain PCA from a gene list
Run.PCA <- function(scTypeEval,
                    gene.list = NULL,
                    black.list = NULL,
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
                         ident.name <- names(mat@ident)
                         
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
                                            ident = setNames(list(mat@ident), ident.name),
                                            sample = mat@sample,
                                            key = "PCA")
                         return(rr)
                      })
   
   names(pca.list) <- names(scTypeEval@data)
   
   # add to scTypeEval object
   for(n in names(pca.list)){
      scTypeEval@reductions[[n]] <- pca.list[[n]]
   }
   
   return(scTypeEval)
   
}

#' @title Add a Custom Dimensionality Reduction to an scTypeEval Object
#'
#' @description
#' This function allows the user to insert pre-computed dimensionality reduction (e.g., PCA, UMAP, 
#' t-SNE, or any embedding) into an \code{scTypeEval} object. The embeddings are stored in the 
#' \code{reductions} slot as a \code{DimRed} object, enabling integration with the scTypeEval 
#' workflow for downstream analysis.
#'
#' @param scTypeEval A \code{scTypeEval} object where the dimensionality reduction will be stored.
#' @param embeddings A numeric matrix of dimension-reduced embeddings (cells/samples x components).
#' @param aggregation Character. Aggregation level of the embeddings. Options: 
#' \code{"single-cell"} (no aggregation) or \code{"pseudobulk"} (aggregated per sample and cell type).
#' @param ident Required. A vector of cell identities (e.g., cell type annotation). 
#' If \code{NULL}, the function attempts to infer it.
#' @param ident.name Character. Name assigned to the provided \code{ident} grouping (default: \code{"custom"}).
#' @param sample Required. A vector indicating sample identity of each observation. 
#' If \code{NULL}, the function attempts to infer it.
#' @param key Optional. Character. Key or label assigned to this dimensionality reduction (e.g., \code{"PCA"}, \code{"UMAP"}).
#' @param gene.list Optional. Character vector or named list of genes associated with the embeddings 
#' (e.g., input features used for the dimensionality reduction). Default: \code{NULL}.
#' @param black.list Optional. Character vector of genes excluded from the dimensionality reduction.
#' Default: \code{NULL}.
#' @param feature.loadings Optional. Matrix of feature loadings corresponding to the embeddings 
#' (e.g., PCA rotation matrix). Default: \code{NULL}.
#' @param filter Logical. If \code{TRUE}, filters cells before adding the embeddings based on 
#' \code{min.samples} and \code{min.cells}. Only supported for \code{"single-cell"} aggregation. 
#' Default: \code{FALSE}.
#' @param min.samples Integer. Minimum number of samples required for retaining a cell type in 
#' single-cell filtering. Default: \code{5}.
#' @param min.cells Integer. Minimum number of cells required per group for filtering. Default: \code{10}.
#' @param verbose Logical. Whether to print messages during execution. Default: \code{TRUE}.
#'
#' @return The modified \code{scTypeEval} object with the new dimensionality reduction stored 
#' in the \code{reductions} slot.
#'
#' @examples
#' \dontrun{
#' # Add precomputed PCA embeddings into scTypeEval
#' sceval <- add.DimReduction(
#'    scTypeEval = sceval,
#'    embeddings = my_pca_matrix,
#'    aggregation = "single-cell",
#'    ident = celltype_annotation,
#'    sample = sample_id,
#'    key = "PCA"
#' )
#' }
#'
#' @seealso \link{Run.PCA}
#'
#' @export add.DimReduction


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
                      ident = setNames(list(ident), ident.name),
                      sample = sample,
                      key = key)
   
   # add to scTypeEval object
   scTypeEval@reductions[[aggregation]]<- rr
   
   return(scTypeEval)
   
}

#' @title Run Dissimilarity Analysis
#'
#' @description
#' Computes dissimilarity between cell populations or pseudobulk profiles 
#' stored in a \code{scTypeEval} object, using one of the supported methods. 
#' Dissimilarity can be computed either on dimensional reduction embeddings 
#' (if available) or on processed gene expression data.
#'
#' @param scTypeEval A \code{scTypeEval} object containing processed data 
#'        and/or dimensional reduction assays.
#' @param method Character. Dissimilarity method to use. Must be one of 
#'        \code{names(dissimilarity_methods)}. Default is \code{"WasserStein"}.
#' @param reduction Logical. Whether to compute dissimilarity on dimensional 
#'        reduction embeddings (if available). Default is \code{TRUE}.
#' @param gene.list Optional. Character vector of genes to include. If 
#'        \code{NULL}, the method will use the default or inherited gene list.
#' @param black.list Optional. Character vector of genes to exclude. If 
#'        \code{NULL}, the method will use the default or inherited blacklist.
#' @param ReciprocalClassifier Character. Classifier to use for RecipClassif 
#'        dissimilarity methods. Default is 
#'        \code{"SingleR"}.
#' @param ncores Integer. Number of cores for parallelization. Default is \code{1}.
#' @param bparam Optional. BiocParallel parameter object to control 
#'        parallelization. Default is \code{NULL}.
#' @param progressbar Logical. Whether to display a progress bar during 
#'        computation. Default is \code{FALSE}.
#' @param verbose Logical. Whether to print progress messages. Default is \code{TRUE}.
#'
#' @details
#' The function supports multiple dissimilarity strategies:
#' \itemize{
#'   \item \code{"Pseudobulk:<distance>"} – computes pairwise distances between 
#'         pseudobulk profiles using the specified distance metric. Supported distances are euclidean, cosine, and pearson.
#'   \item \code{"WasserStein"} – computes Wasserstein distances between groups 
#'         of embeddings or cells.
#'   \item \code{"RecipClassif:<method>"} – assigns cells pairwise between samples using 
#'         the specified classifier, then computes dissimilarity between assignments.
#'         Supported methods are 'match' and 'score'.
#' }
#'
#' If \code{reduction = TRUE}, the function expects that dimensional reduction 
#' embeddings have been added previously via \code{Run.PCA()} or 
#' \code{add.DimReduction()}. If unavailable, set \code{reduction = FALSE} 
#' to compute dissimilarity on processed expression data instead.
#'
#' @return
#' An updated \code{scTypeEval} object with a new \code{DissimilarityAssay} 
#' stored in \code{scTypeEval@dissimilarity[[method]]}.
#'
#' @examples
#' \dontrun{
#' # Run Wasserstein dissimilarity on dimensional reduction embeddings
#' scTypeEval <- Run.Dissimilarity(scTypeEval,
#'                                method = "WasserStein",
#'                                reduction = TRUE)
#'
#' # Run pseudobulk Euclidean dissimilarity on processed data
#' scTypeEval <- Run.Dissimilarity(scTypeEval, 
#'                                 method = "Pseudobulk:euclidean", 
#'                                 reduction = FALSE)
#'                                 
#' # Run RecipClassif classification-based dissimilarity with SingleR
#' scTypeEval <- Run.Dissimilarity(scTypeEval, 
#'                                 method = "RecipClassif:Match", 
#'                                 ReciprocalClassifier = "SingleR")
#' }
#'
#' @seealso 
#' \code{\link{add.DimReduction}}, \code{\link{Run.PCA}}, 
#' \code{\link{Run.ProcessingData}}
#'
#' @export Run.Dissimilarity



Run.Dissimilarity <- function(scTypeEval,
                              method = "WasserStein",
                              reduction = TRUE,
                              gene.list = NULL,
                              black.list = NULL,
                              ReciprocalClassifier = "SingleR",
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
   
   if(reduction && method %in% no_dr_ds){
      warning("No dimensional reduction dissimilarity computation supported for: ",
              paste(no_dr_ds, collapse = ", "), ". Switching to `reduction=FALSE`\n")
      reduction <- FALSE
   }
   
   if(reduction){
      mat_ident <- scTypeEval@reductions[[slot]]
      if(is.null(mat_ident)){
         stop("No dimensional reduction slot found for ", slot ,
              ". Please run before `Run.PCA()` or add a valid dimensional reduction assay.\n")
      }
      mat <- mat_ident@embeddings
      ident.name <- names(mat_ident@ident)
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
      ident.name <- names(mat_ident@ident)
      mat_ident <- .general_filtering(mat_ident,
                                      black.list = black.list,
                                      gene.list = gene.list,
                                      verbose = verbose)
      mat <- mat_ident@matrix
   }
   
   # inheret ident and sample
   if(slot == "single-cell") {
      group_levels <- levels(mat_ident@group)
      ident <- sapply(group_levels, function(x){strsplit(x, "_")[[1]][2]}) |>
         factor()
      ident <- setNames(list(ident), ident.name)
      sample <- sapply(group_levels, function(x){strsplit(x, "_")[[1]][1]}) |>
         factor()
      
   } else if(slot == "pseudobulk"){
      ident <- mat_ident@ident
      if(!is.list(ident)){
         ident <- setNames(list(ident), ident.name)
      }
      sample <- mat_ident@sample
   }
   
   
   # set paralelization
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   # split matrix per sample for ReciprocalClassif, and for celltype & sample for Wasserstein, no need for the rest
   
   aggregation <- strsplit(method, ":")[[1]][1] |> tolower()
   dist.type <- strsplit(method, ":")[[1]][2] |> tolower()
   
   dist <- switch(
      aggregation,
      "pseudobulk" = get.distance(norm.mat = mat,
                                  distance.method = dist.type,
                                  verbose = verbose),
      "wasserstein" = compute_wasserstein(mat = mat,
                                          group = mat_ident@group,
                                          bparam = param,
                                          verbose = verbose),
      "recipclassif" = RecipClassif(mat = mat,
                                    ident = ident[[1]],
                                    sample = sample,
                                    group = mat_ident@group,
                                    classifier = ReciprocalClassifier,
                                    method = dist.type,
                                    bparam = param,
                                    verbose = verbose),
      stop(aggregation, "is not a supported method.\n")
   )
   
   # build dissimilarity object
   rr <- methods::new("DissimilarityAssay",
                      dissimilarity = dist,
                      method = method,
                      gene.list = gene.list,
                      black.list = black.list,
                      aggregation = aggregation,
                      ident = ident,
                      sample = sample)
   
   scTypeEval@dissimilarity[[method]] <- rr
   
   
   return(scTypeEval)
   
}



#' @title Compute Consistency of Cell Type Annotations Across samples
#'
#' @description
#' Computes internal validation metrics (consistency measures) for cell type 
#' annotations based on dissimilarity assays stored in a \code{scTypeEval} object. 
#' For each dissimilarity representation, one or more internal validation metrics 
#' are calculated per cell type and returned in a tidy \code{data.frame}.
#'
#' @param scTypeEval A \code{scTypeEval} object containing one or more 
#'        dissimilarity assays (see \code{Run.Dissimilarity}).
#' @param dissimilarity.slot Character. Which dissimilarity assay(s) to use. 
#'        Can be \code{"all"} (default) or the name(s) of specific slots stored 
#'        in \code{scTypeEval@dissimilarity}.
#' @param Consistency.metric Character vector. Internal validation metrics to compute. 
#'        Supported options include:
#'        \itemize{
#'          \item \code{"silhouette"} – cohesion/separation score
#'          \item \code{"2label.silhouette"} – Variant of the silhouette score where,
#'           for each cell, the within-cluster distance is compared against the average
#'           distance to all cells outside its cluster (rather than just the closest
#'           other cluster). This effectively treats the clustering as a two-group
#'           problem: "own cell type" vs. "all others." Scores are then averaged per
#'           cell type. Higher values indicate stronger separation of a cell type
#'           from the rest of the dataset as a whole.
#'          \item \code{"NeighborhoodPurity"} – For each cell, computes the fraction of its K nearest
#'          neighbors that belong to the same cell type. Scores are then averaged per cell type.
#'          Optionally, values can be normalized relative to the expected proportion by chance.
#'          Higher scores indicate that cells are surrounded by neighbors of the same type,
#'          reflecting strong local label consistency.
#'          \item \code{"ward.PropMatch"} – For each true cell type, identifies the cluster that contains
#'           the largest number of its cells (the dominant cluster) and computes the proportion of cells
#'           from that cell type that fall into this cluster. Optionally, this proportion can be normalized
#'           relative to the expected proportion by chance. Higher values indicate better alignment between
#'           true labels and cluster assignments.
#'          \item \code{"Orbital.medoid"} – For each cell type, identifies a representative medoid cell
#'           (the cell minimizing the total distance to all other cells of the same type).
#'            Then, for each non-medoid cell, it checks whether the cell is closer to its own medoid
#'             than to the medoids of other cell types. The metric for each cell type is the proportion
#'             of cells that are closer to their own medoid than to any other medoid. 
#'              Optionally, this proportion can be normalized relative to the expected proportion by chance.
#'              Higher values indicate that cells are well-clustered around their medoid.
#'          \item \code{"Average.similarity"} – Measures how similar cells are within
#'           the same cell type relative to cells of other types. For each cell,
#'           it computes the average distance to other cells in its group and to
#'          cells outside its group, then combines these into a normalized score (higher = better).
#'          Scores are then averaged per cell type. This metric is similar in spirit to a one-label silhouette.
#'        }
#'        Default: all supported metrics.
#' @param KNNGraph_k Integer. Number of nearest neighbors to use for 
#'        graph-based metrics (e.g. \code{"NeighborhoodPurity"}). Default is \code{5}.
#' @param hclust.method Character. Agglomeration method passed to \code{\link{hclust}} 
#'        for Ward-based metrics. Default is \code{"ward.D2"}.
#' @param normalize Logical. Whether to normalize metric values (e.g., scaling 
#'        across dissimilarity methods). Default is \code{FALSE}.
#' @param return.scTypeEval Logical. Whether to return data frame with inter-sample consistencies or store within scTypeEval@consistency slot. Default is \code{FALSE}.
#' @param verbose Logical. Whether to print progress messages. Default is \code{TRUE}.
#'
#' @details
#' This function builds upon the dissimilarity assays generated by 
#' \code{\link{Run.Dissimilarity}}. For each selected dissimilarity representation, 
#' the chosen internal validation metrics are computed and stored in a long-format 
#' data frame, allowing downstream comparison across cell types, metrics, 
#' and dissimilarity methods.
#'
#' @return 
#' A \code{data.frame} with the following columns:
#' \itemize{
#'   \item \code{celltype} – the annotation/group label
#'   \item \code{measure} – numeric consistency score
#'   \item \code{consistency.metric} – the metric name
#'   \item \code{dissimilarity_method} – the dissimilarity method used
#'   \item \code{ident} – the identity class (from \code{scTypeEval@ident})
#' }
#'
#' @examples
#' \dontrun{
#' # Compute silhouette and Neighborhood Purity consistency metrics
#' consistency_df <- get.Consistency(scTypeEval,
#'                                   dissimilarity.slot = "WasserStein",
#'                                   Consistency.metric = c("silhouette", "NeighborhoodPurity"))
#'
#' # Run all metrics on all available dissimilarity slots
#' consistency_df <- get.Consistency(scTypeEval, dissimilarity.slot = "all")
#' }
#'
#' @seealso 
#' \code{\link{Run.Dissimilarity}}
#'
#' @export get.Consistency


get.Consistency <- function(scTypeEval,
                            dissimilarity.slot = "all", 
                            Consistency.metric = c("silhouette",
                                                   "2label.silhouette",
                                                   "NeighborhoodPurity",
                                                   "ward.PropMatch",
                                                   "Orbital.medoid",
                                                   "Average.similarity"),
                            KNNGraph_k = 5,
                            hclust.method = "ward.D2",
                            normalize = FALSE,
                            return.scTypeEval = FALSE,
                            verbose = TRUE
                            
){
   
   diss.assays <- .check_dissimilarityAssays(scTypeEval, slot = dissimilarity.slot)
   
   consist.list <- lapply(diss.assays,
                          function(da){
                             assay <- scTypeEval@dissimilarity[[da]]
                             
                             dist <- assay@dissimilarity
                             # if RecipClassif match, convert 0.5 (for plotting to 1)
                             if (da == "RecipClassif:Match") {
                                dist[dist == 0.5] <- 1
                             }
                             
                             ident <- assay@ident[[1]]
                             ident.name <- names(assay@ident)
                             # all expected clusters
                             all_clusters <- purge_label(unique(scTypeEval@metadata[[ident.name]]))
                             all_clusters <- all_clusters[!is.na(all_clusters)]
                             
                             # compute internal validation metrics
                             if(verbose){message("Computing internal validation metrics for ", da, " ... \n")}
                             con <- calculate_IntVal_metric(dist = dist,
                                                            metrics = Consistency.metric,
                                                            ident = ident,
                                                            KNNGraph_k = KNNGraph_k,
                                                            hclust.method = hclust.method,
                                                            normalize = normalize)
                             dfl <- lapply(names(con),
                                           function(int){
                                              # ensure all clusters/celltypes are present
                                              # cell type no passing the threshold will have a value of 0
                                              con[[int]] <- con[[int]][all_clusters]
                                              con[[int]][is.na(con[[int]])] <- 0
                                              names(con[[int]]) <- all_clusters
                                              # start consistency dataframe object
                                              r <- data.frame(celltype = names(con[[int]]),
                                                              measure = con[[int]],
                                                              consistency.metric = int)
                                              return(r)
                                           })
                             
                             cons <- do.call(rbind, dfl) |>
                                dplyr::mutate(dissimilarity_method = assay@method,
                                              ident = names(assay@ident))
                             
                             return(cons)
                          })
   
   consist <- do.call(rbind, consist.list) 
   
   if(!return.scTypeEval){
      return(consist)
   } else {
      ident.name <- unique(consist$ident)
      scTypeEval@consistency[[ident.name]] <- consist
      return(scTypeEval)
   }
   
}


#' @title Plot PCA Results from scTypeEval Object
#' 
#' @description
#' This function visualizes Principal Component Analysis (PCA) results stored in the 
#' \code{reductions} slot of an \code{scTypeEval} object.
#'
#' @param scTypeEval An \code{scTypeEval} object containing PCA results in the \code{reductions} slot.
#' @param reduction.slot Character. Name(s) of the reduction(s) to plot. If \code{"all"} 
#' (default), all available PCA reductions in the object are plotted.
#' @param label Logical. Whether to add cluster labels to the PCA plot (default: \code{TRUE}).
#' @param dims Integer vector of length 2. The principal component (PC) dimensions to plot 
#' (default: \code{c(1,2)}).
#' @param show.legend Logical. Whether to display a legend (default: \code{FALSE}).
#'
#' @return A named list of PCA plots (\link[ggplot2]{ggplot} objects) corresponding to 
#' the PCA analyses stored in the \code{reductions} slot of an \code{scTypeEval} object 
#' (generated by \link{add.PCA}). If only one reduction is selected, a single ggplot object 
#' is returned.
#'
#' @seealso \link{add.PCA}
#'
#' @examples
#' \dontrun{
#' # plot all PCA dimensional reductions
#' pca_plots <- plot.PCA(sceval) # scTypeEval object with previously run `add.PCA()`
#' }
#'
#' @export plot.PCA



plot.PCA <- function(scTypeEval,
                     reduction.slot = "all",
                     label = TRUE,
                     dims = c(1,2),
                     show.legend = FALSE) {
   
   red.assays <- .check_DimRedAssays(scTypeEval, slot = reduction.slot)
   
   pls <- lapply(red.assays,
                 function(da){
                    assay <- scTypeEval@reductions[[da]]
                    ident <- assay@ident[[1]]
                    title <- paste(da, names(assay@ident), sep = " - ")
                    
                    df <- t(assay@embeddings[dims,]) |>
                       as.data.frame() |>
                       dplyr::mutate(ident = ident)
                    
                    # compute variance of PCs
                    vrs <- var_PCA(t(assay@embeddings))[dims]
                    vrs <- round(vrs*100, 2)
                    
                    labs <- paste0("PC", dims, " (", vrs, "%)")
                    
                    pl <- helper.plot.scatter(df,
                                              show.legend = show.legend,
                                              label = label) +
                       ggplot2::labs(x = labs[1],
                                     y = labs[2],
                                     title = title)
                    return(pl)
                 })
   
   names(pls) <- red.assays
   
   if(length(pls) == 1){
      pls <- unlist(pls)
   }
   
   return(pls)
}

#' @title Perform Hierarchical Clustering on scTypeEval Dissimilarity Matrices
#'
#' @description
#' This function performs hierarchical clustering using precomputed dissimilarity
#' matrices stored in the \code{dissimilarity} slot of an \code{scTypeEval} object.
#' Each dissimilarity assay is clustered independently using the specified hierarchical
#' clustering method.
#'
#' @param scTypeEval An \code{scTypeEval} object containing dissimilarity matrices in the \code{dissimilarity} slot.
#' @param dissimilarity.slot Character string. Specifies which dissimilarity assays to cluster.
#'   Use \code{"all"} (default) to include all available assays.
#' @param hierarchy.method Character string specifying the hierarchical clustering method
#'   (default: \code{"ward.D2"}). See \link[stats]{hclust} for available options.
#' @param verbose Logical. If \code{TRUE}, prints messages about progress (default: \code{TRUE}).
#'
#' @return A list of clustering results, one per dissimilarity assay. Each element
#'   contains a contingency table of cluster assignments versus input labels.
#'
#' @details
#' For each dissimilarity assay, the function applies \link[stats]{hclust}
#' to the stored dissimilarity matrix. The number of clusters is set equal
#' to the number of unique identities provided in the assay metadata.
#'
#' @examples
#' \dontrun{
#' # Perform hierarchical clustering on all available dissimilarity matrices
#' hier_results <- get.hierarchy(sceval)
#' }
#'
#' @seealso \link{Run.dissimilarity}, \link[stats]{hclust}
#'
#' @export get.hierarchy


get.hierarchy <- function(scTypeEval,
                          dissimilarity.slot = "all", 
                          hierarchy.method = "ward.D2",
                          verbose = TRUE
){
   diss.assays <- .check_dissimilarityAssays(scTypeEval, slot = dissimilarity.slot)
   
   hier.list <- lapply(diss.assays,
                       function(da){
                          assay <- scTypeEval@dissimilarity[[da]]
                          
                          dist <- assay@dissimilarity
                          ident <- assay@ident[[1]]
                          
                          # compute internal validation metrics
                          if(verbose){message("Computing hierarchical clustering for ", da, " ... \n")}
                          attr(dist, "Labels") <- ident
                          
                          hclust_result <- stats::hclust(dist,
                                                         method = hierarchy.method)
                          clusters <- stats::cutree(hclust_result,
                                                    k = length(unique(ident)))
                          
                          t <- table(clusters, names(clusters))
                          return(t)
                       })
   
   names(hier.list) <- diss.assays
   
   if(length(hier.list) == 1){
      hier.list <- unlist(hier.list)
   }
   
   return(hier.list)
}

#' @title Compute K-Nearest Neighbor (KNN) Composition from Dissimilarity Assays
#'
#' @description
#' This function computes KNN-based neighborhood composition scores from dissimilarity 
#' matrices stored in a `scTypeEval` object. For each dissimilarity assay, it constructs 
#' a KNN graph, computes the cell type composition among neighbors 
#' and aggregates results at the group level.
#'
#' @param scTypeEval An object containing dissimilarity matrices and metadata.
#' @param dissimilarity.slot A character string specifying which dissimilarity assay(s) 
#'   to use. If `"all"`, all available dissimilarity assays are processed (default: `"all"`).
#' @param KNNGraph_k Integer; the number of neighbors to consider in the KNN graph (default: `5`).
#' @param normalize Logical; if `TRUE`, normalizes neighbor proportions relative to expected 
#'   frequencies of each cell type (default: `FALSE`).
#' @param verbose Logical; if `TRUE`, prints progress messages during computation (default: `TRUE`).
#'
#' @return 
#' A list (or a single data frame if only one assay is processed) where each element 
#' corresponds to a dissimilarity assay. Each element contains a data frame with rows 
#' corresponding to reference cell types and columns representing the mean proportion of 
#' neighbors belonging to each cell type.
#'
#' @details
#' The function extracts dissimilarity matrices from the `scTypeEval` object and applies 
#' a KNN graph construction. For each assay, it computes neighbor 
#' cell-type proportions per cell and then aggregate them by cell type. If 
#' normalization is enabled, the observed neighbor proportions are scaled relative to 
#' the expected global frequency of each cell type.
#'
#' @examples
#' \dontrun{
#' result <- get.NN(scTypeEval = sceval,
#'                  KNNGraph_k = 5,
#'                  normalize = TRUE)
#' }
#'
#' @export get.NN


get.NN <- function(scTypeEval,
                   dissimilarity.slot = "all",
                   KNNGraph_k = 5,
                   normalize = FALSE,
                   verbose = TRUE
){
   diss.assays <- .check_dissimilarityAssays(scTypeEval, slot = dissimilarity.slot)
   
   nn.list <- lapply(diss.assays,
                     function(da){
                        assay <- scTypeEval@dissimilarity[[da]]
                        
                        dist <- assay@dissimilarity
                        ident <- assay@ident[[1]]
                        
                        # compute internal validation metrics
                        if(verbose){message("Computing hierarchical clustering for ", da, " ... \n")}
                        attr(dist, "Labels") <- ident
                        
                        nn <- nn.helper(dist = dist,
                                        ident = ident,
                                        KNNGraph_k = KNNGraph_k,
                                        normalize = normalize)
                        
                        return(nn)
                     })
   
   names(nn.list) <- diss.assays
   
   if(length(nn.list) == 1){
      nn.list <- unlist(nn.list)
   }
   
   return(nn.list)
}

#' @title Plot MDS Results from scTypeEval Object
#'
#' @description
#' This function visualizes Multidimensional Scaling (MDS) results computed from 
#' dissimilarity assays stored in an `scTypeEval` object. For each dissimilarity assay, 
#' MDS is performed using `stats::cmdscale()` and results are displayed as scatterplots.
#'
#' @param scTypeEval An `scTypeEval` object containing one or more dissimilarity assays.
#' @param dissimilarity.slot Character string specifying which dissimilarity assay(s) 
#'   to use. If `"all"`, all available dissimilarity assays are plotted (default: `"all"`).
#' @param label Logical; whether to add cluster labels to the MDS plot (default: `TRUE`).
#' @param dims Integer vector of length 2; the MDS dimensions to plot (default: `c(1, 2)`).
#' @param show.legend Logical; whether to display a legend (default: `FALSE`).
#'
#' @return 
#' A named list of MDS plots (`ggplot2` objects), one per dissimilarity assay. 
#' If only a single assay is processed, a single `ggplot2` object is returned.
#'
#' @details
#' For each selected dissimilarity assay, the function:
#' \enumerate{
#'   \item Extracts the dissimilarity matrix and cell-type identities.
#'   \item Performs classical multidimensional scaling (MDS) with `cmdscale()`.
#'   \item Generates a 2D scatterplot with group labels and optional legends using 
#'   an internal plotting helper.
#' }
#' The plot axes correspond to the requested MDS dimensions (`dims`).
#'
#' @examples
#' \dontrun{
#' # Plot MDS for all dissimilarity assays
#' mds_plots <- plot.MDS(sceval)
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
                    ident <- assay@ident[[1]]
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
      pls <- pls[[1]]
   }
   
   return(pls)
}

#' @title Plot Heatmaps of Dissimilarity Matrices
#'
#' @description
#' Visualize dissimilarity matrices stored in an `scTypeEval` object as annotated
#' heatmaps. Each selected dissimilarity assay is shown as a heatmap where rows and
#' columns represenent cell type and sample, ordered by cell type and optionally sorted by similarity
#' or consistency metrics. Group boundaries are marked to highlight cell-type consistency.
#'
#' @param scTypeEval An `scTypeEval` object containing one or more dissimilarity assays.
#' @param dissimilarity.slot Character string specifying which dissimilarity assay(s) 
#'   to plot. If `"all"`, all available dissimilarity assays are included (default: `"all"`).
#' @param sort.similarity Optional. Character string naming a dissimilarity assay to use
#'   for ordering cells by similarity (hierarchical clustering within each cell type).
#' @param sort.consistency Optional. Character string specifying a consistency metric
#'   (passed to `get.Consistency()`) for ordering cell types by overall consistency.
#' @param low.color Color for the low end of the heatmap gradient (default: `"black"`).
#' @param high.color Color for the low end of the heatmap gradient (default: `"white"`).
#' @param hclust.method Clustering method to use when `sort.similarity` is provided
#'   (default: `"ward.D2"`).
#' @param verbose Logical. Whether to print progress and diagnostic messages (default: `TRUE`).
#' @param ... Additional arguments passed to `get.Consistency()`.
#'
#' @return
#' A named list of `ggplot2` objects, one per dissimilarity assay.
#' If only a single assay is selected, the corresponding heatmap plot is returned directly.
#'
#' @details
#' Ordering logic:
#' \itemize{
#'   \item If both `sort.similarity` and `sort.consistency` are `NULL`, cell types are
#'         ordered alphabetically, and cells within each type are ordered alphabetically.
#'   \item If only `sort.consistency` is provided, cell types are ordered by the
#'         selected consistency metric, and cells within each type are ordered alphabetically.
#'   \item If only `sort.similarity` is provided, cell types are ordered by hierarchical
#'         clustering of average similarities, and cells within each type are clustered.
#'   \item If both are provided, cell types are ordered by consistency, while cells
#'         within each type are ordered by similarity clustering.
#' }
#'
#' Each heatmap is rendered with ggplot2, with group boundaries and axis labels
#' indicating cell-type structure.
#'
#' @examples
#' \dontrun{
#' # Plot heatmaps for all dissimilarity assays
#' hm_all <- plot.Heatmap(sceval)
#'
#' # Plot with cell types ordered by consistency
#' hm_cons <- plot.Heatmap(sceval, sort.consistency = "silhouette")
#'
#' # Plot with cells clustered by similarity within cell type
#' hm_sim <- plot.Heatmap(sceval, sort.similarity = "WasserStein")
#'
#' # Plot with consistency ordering of cell types and similarity ordering within
#' hm_both <- plot.Heatmap(sceval, sort.consistency = "silhouette", sort.similarity = "WasserStein")
#' }
#' 
#' @seealso \link{Run.dissimilarity}, \link[stats]{hclust}, \link{get.Consistency}
#'
#' @export plot.Heatmap



plot.Heatmap <- function(scTypeEval,
                         dissimilarity.slot = "all",
                         sort.similarity = NULL,
                         sort.consistency = NULL,
                         low.color = "black",
                         high.color = "white",
                         hclust.method = "ward.D2",
                         verbose = TRUE,
                         ...
) {
   
   # check if dissimilarity object are present in object
   diss.assays <- .check_dissimilarityAssays(scTypeEval, slot = dissimilarity.slot)
   
   # if not indicated just order cell types based on alphabetical
   if(is.null(sort.similarity) && is.null(sort.consistency)){
      if(verbose){message("No ordering (sort.similarity or sort.consistency) indicated,
                          sorting cell type by alphabetical order.")}
   }
   
   # Compute similarity of cell types if indicated
   if(!is.null(sort.similarity)){
      # get dissimilarity matrix
      diss.assays.cluster <- .check_dissimilarityAssays(scTypeEval, slot = sort.similarity)
      assay_cluster <- scTypeEval@dissimilarity[[diss.assays.cluster]]
      d_mat_cluster <- as.matrix(assay_cluster@dissimilarity)
      # get similarities between cell types
      if(verbose){message("Computing cell type similarity based on ", diss.assays.cluster, ".")}
      ident_levels <- levels(assay_cluster@ident[[1]])
      mat <- matrix(0, nrow = length(ident_levels),
                    ncol = length(ident_levels),
                    dimnames = list(ident_levels, ident_levels))
      for (i in ident_levels) {
         for (j in ident_levels) {
            sub <- d_mat_cluster[assay_cluster@ident[[1]] == i,
                                 assay_cluster@ident[[1]] == j,
                                 drop = FALSE]
            mat[i, j] <- mean(sub)
         }
      }
      hc <- stats::hclust(as.dist(mat),
                          method = hclust.method)
      ident_order <- hc$labels[hc$order]
   } else {
      ident_order <- NULL
   }
   
   pls <- lapply(diss.assays,
                 function(da){
                    assay <- scTypeEval@dissimilarity[[da]]
                    title <- paste(da, names(assay@ident), sep = " - ")
                    
                    d_mat_plot <- as.matrix(assay@dissimilarity)
                    
                    annot_ident <- assay@ident[[1]]
                    
                    #  order by consistency if indicated
                    if (!is.null(sort.consistency)) {
                       
                       if(verbose){message("Computing consistency metric for ", sort.consistency, ".\n")}
                       
                       consis <- get.Consistency(scTypeEval,
                                                 dissimilarity.slot = da,
                                                 Consistency.metric = sort.consistency,
                                                 verbose = F,
                                                 ...) |>
                          dplyr::arrange(measure) |>
                          dplyr::pull(measure, name = celltype)
                       
                       if(verbose){print(consis)}
                       
                       ident_order <- names(consis)  
                    } 
                    
                    if(is.null(ident_order)){
                       # order alphabetically
                       ident_order <- sort(levels(annot_ident))
                    }
                    
                    # reorder cells within cell type: first by ident order, then cluster within ident
                    if(!is.null(sort.similarity)){
                       ordering <- unlist(lapply(ident_order, function(iden) {
                          cells <- names(annot_ident)[annot_ident == iden]
                          if (length(cells) > 2) {
                             hc <- stats::hclust(as.dist(d_mat_cluster[cells, cells]),
                                                 method = hclust.method)
                             cells[hc$order]
                          } else {
                             cells
                          }
                       }))
                    } else {
                       ordering <- unlist(lapply(ident_order, function(iden) {
                          cells <- names(annot_ident)[annot_ident == iden] |>
                             sort()
                          return(cells)
                       }))
                    }
                    
                    # reorder matrix + annotation
                    d_mat_plot <- d_mat_plot[ordering, ordering, drop = F]
                    
                    # convert matrix to tidy df
                    df <- as.data.frame(as.table(d_mat_plot))
                    colnames(df) <- c("row", "col", "value")
                    
                    # ident boundaries + label positions
                    counts <- as.numeric(table(annot_ident)[ident_order])
                    sep_lines <- cumsum(counts)
                    ident_labels <- ident_order
                    
                    # midpoints for axis ticks
                    x_ticks <- sep_lines - counts/2 + 0.5
                    y_ticks <- length(ordering) - (sep_lines - counts/2 + 0.5) + 1
                    
                    # build plot
                    p <- ggplot2::ggplot(df, ggplot2::aes(x = col, y = row, fill = value)) +
                       ggplot2::geom_tile() +
                       ggplot2::scale_fill_gradient(low = low.color,
                                                    high = high.color) +
                       ggplot2::scale_x_discrete(
                          breaks = ordering[x_ticks],
                          labels = ident_labels,
                          expand = c(0, 0)
                       ) +
                       ggplot2::scale_y_discrete(
                          breaks = rev(ordering)[y_ticks],
                          labels = ident_labels,
                          expand = c(0, 0)
                       ) +
                       ggplot2::theme_minimal(base_size = 12) +
                       ggplot2::theme(
                          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
                          axis.text.y = ggplot2::element_text(hjust = 1),
                          axis.ticks = ggplot2::element_blank(),
                          axis.title = ggplot2::element_blank(),
                          legend.position = "none",
                          panel.grid = ggplot2::element_blank()
                       ) +
                       ggplot2::geom_hline(yintercept = sep_lines + 0.5, color = "grey70", size = 0.4) +
                       ggplot2::geom_vline(xintercept = sep_lines + 0.5, color = "grey70", size = 0.4) +
                       ggplot2::labs(title = title)
                    
                    return(p)
                 })
   
   names(pls) <- diss.assays
   
   if(length(pls) == 1){
      pls <- pls[[1]]
   }
   
   return(pls)
}


#' @title Load a Single-Cell Object from File
#'
#' @description
#' Loads single-cell datasets from `.rds` or `.h5ad` files into R.
#' Supports \pkg{Seurat}, \pkg{SingleCellExperiment}, and \pkg{anndata} objects.
#' Depending on the input type and the `split` parameter, the function
#' either returns the original object or extracts and returns the
#' raw counts matrix and metadata.
#'
#' @param path Character string. Path to the file to load.
#'   Supported formats: `.rds` (Seurat or SingleCellExperiment objects)
#'   and `.h5ad` (AnnData objects).
#' @param split Logical (default: `TRUE`).  
#'   If `TRUE`, returns a list with two elements:
#'   \itemize{
#'     \item \code{counts}: A sparse \code{dgCMatrix} of raw counts.
#'     \item \code{metadata}: A \code{data.frame} of cell metadata.
#'   }
#'   If `FALSE`, returns the loaded object as-is
#'   (Seurat, SingleCellExperiment, or AnnData).
#'
#' @return
#' If `split = TRUE`, a list containing:
#' \itemize{
#'   \item \code{counts} – sparse counts matrix
#'   \item \code{metadata} – cell metadata
#' }  
#' If `split = FALSE`, returns the loaded single-cell object directly.
#'
#' @details
#' - **`.rds` input**:  
#'   - If the object is a \pkg{Seurat} object, the raw counts are extracted
#'     from the "RNA" assay, and cell metadata is taken from `object@meta.data`.
#'   - If the object is a \pkg{SingleCellExperiment}, counts are extracted from
#'     the `"counts"` assay, and cell metadata is taken from `colData(object)`.
#'   - Other `.rds` object types are not supported.
#'
#' - **`.h5ad` input**:  
#'   Requires the \pkg{anndata} package. The `"counts"` layer must be present,
#'   otherwise the function will stop with an error. Counts are transposed
#'   to match R’s cell-by-gene convention.
#'
#' @examples
#' \dontrun{
#' # Load as counts and metadata
#' sce_split <- load_singleCell_object("dataset.rds", split = TRUE)
#'
#' # Load Seurat object as-is
#' seurat_obj <- load_singleCell_object("dataset.rds", split = FALSE)
#'
#' # Load AnnData (.h5ad) file and extract counts + metadata
#' adata_split <- load_singleCell_object("dataset.h5ad", split = TRUE)
#' }
#'
#' @export load_singleCell_object


load_singleCell_object <- function(path,
                                   split = TRUE) {
   
   if (!file.exists(path)) stop("File does not exist: ", path)
   
   # Initialize
   object <- NULL
   counts <- NULL
   metadata <- NULL
   
   if (grepl("rds$", path, ignore.case = T)) {
      object <- readRDS(path)
      
      if (inherits(object, "Seurat")) {
         if (!requireNamespace("Seurat", quietly = TRUE)) {
            stop("The 'Seurat' package is required to read Seurat objects. Please install it.")
         }
         rcounts <- Seurat::GetAssayData(object,
                                         assay = "RNA",
                                         layer = "counts")
         counts <- as(rcounts, "dgCMatrix")
         metadata <- as.data.frame(object@meta.data)
         
      } else if (inherits(object, "SingleCellExperiment")) {
         if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
            stop("The 'SummarizedExperiment' package is required to read SCE objects. Please install it.")
         }
         rcounts <- SummarizedExperiment::assay(object, "counts")
         counts <- as(rcounts, "dgCMatrix")
         metadata <- as.data.frame(SummarizedExperiment::colData(object))
         
      } else {
         stop("Unsupported .rds object type: ", class(object))
      }
      
   } else if (grepl("h5ad$", path, ignore.case = T)) {
      if (!requireNamespace("anndata", quietly = TRUE)) {
         stop("The 'anndata' package is required to read .h5ad files. Please install it.")
      }
      
      object <- anndata::read_h5ad(path)
      
      if ("counts" %in% names(object$layers)) {
         message("Using 'counts' layer from .h5ad.")
         rcounts <- object$layers[["counts"]]
         
      } else if (!is.null(object$X)) {
         message("No 'counts' layer found. Using main matrix 'X'.")
         rcounts <- object$X
         
      } else if (length(object$layers) > 0) {
         first_layer <- names(object$layers)[1]
         message("No 'counts' layer or 'X' found. Using first available layer: ", first_layer)
         rcounts <- object$layers[[first_layer]]
         
      } else {
         stop("No usable expression matrix found in the .h5ad file. 
               Searched for layers[['counts']], X, or any layer.")
      }
      
      rcounts <- Matrix::t(object$layers[["counts"]])
      counts <- as(as(rcounts, "CsparseMatrix"), "dgCMatrix")
      metadata <- as.data.frame(object$obs)
      
   } else {
      stop("Unsupported file format. Please provide a .rds or .h5ad file.")
   }
   
   if (split) {
      return(list(counts = counts,
                  metadata = metadata))
   } else {
      return(object)
   }
}


#' @title Wrapper Function to Compute Dissimilarities from Single-Cell Data
#'
#' @description
#' A high-level convenience function that initializes an \code{scTypeEval} object 
#' from a count matrix and metadata, performs preprocessing (normalization, 
#' filtering, optional dimensionality reduction), defines gene lists, 
#' and computes one or more dissimilarity metrics between cell populations.
#' 
#' This function integrates multiple internal \pkg{scTypeEval} pipeline steps, 
#' including data preparation, HVG selection, PCA reduction, and dissimilarity computation, 
#' providing a streamlined workflow for single-cell data evaluation.
#'
#' @param count_matrix A numeric or sparse \code{dgCMatrix} of raw counts 
#'   (genes as rows, cells as columns). Seurat and SingleCellExperiment object are accepted too.
#' @param metadata A \code{data.frame} containing cell-level metadata.
#' @param ident Character string indicating the metadata column specifying 
#'   cell identities (e.g., cell types or clusters).
#' @param sample Character string specifying the metadata column containing 
#'   sample identifiers (used for pseudobulk aggregation).
#' @param aggregation Method to group cells, either `"single-cell"` or `"pseudobulk"`. Default is `"single-cell"`.
#' @param gene.list Optional named list of gene sets to include in the analysis.  
#'   If \code{NULL}, highly variable genes (HVGs) are automatically computed.
#' @param reduction Logical; if \code{TRUE}, performs PCA dimensionality 
#'   reduction prior to dissimilarity computation (default: \code{TRUE}).
#' @param ndim Integer; number of principal components to retain 
#'   when \code{reduction = TRUE} (default: \code{30}).
#' @param black.list Optional character vector of genes to exclude from analysis.  
#'   If \code{NULL}, no genes are blacklisted.
#' @param normalization.method Character string specifying the normalization 
#'   method to apply. Options include \code{"Log1p"}, \code{"CLR"}, and \code{"pearson"} 
#'   (default: \code{"Log1p"}).
#' @param dissimilarity.method Character vector of dissimilarity metrics to compute.  
#'   Available options include:
#'   \itemize{
#'     \item \code{"WasserStein"}
#'     \item \code{"Pseudobulk:Euclidean"}
#'     \item \code{"Pseudobulk:Cosine"}
#'     \item \code{"Pseudobulk:Pearson"}
#'     \item \code{"RecipClassif:Match"}
#'     \item \code{"RecipClassif:Score"}
#'   }
#'   Multiple methods can be provided for comparative evaluation.
#' @param min.samples Integer; minimum number of samples required for pseudobulk analysis (default: \code{5}).
#' @param min.cells Integer; minimum number of cells required per group for inclusion (default: \code{10}).
#' @param ncores Integer; number of CPU cores for parallel execution (default: \code{1}).
#' @param bparam Optional \code{BiocParallelParam} object for fine-grained parallelization control.  
#'   Overrides \code{ncores} if provided.
#' @param progressbar Logical; if \code{TRUE}, displays a progress bar (default: \code{FALSE}).
#' @param verbose Logical; if \code{TRUE}, prints progress messages (default: \code{TRUE}).
#'
#' @return
#' An updated \code{scTypeEval} object containing:
#' \itemize{
#'   \item Normalized and filtered data
#'   \item HVG gene sets or user-provided gene lists
#'   \item PCA reductions (if enabled)
#'   \item Computed dissimilarity matrices for each selected method
#' }
#'
#' @details
#' This wrapper combines multiple pipeline steps from \pkg{scTypeEval}:
#' \enumerate{
#'   \item \code{create.scTypeEval()} — initializes the object.
#'   \item \code{Run.ProcessingData()} — performs normalization and filtering.
#'   \item \code{Run.HVG()} or \code{add.GeneList()} — defines gene sets.
#'   \item \code{Run.PCA()} — performs PCA if \code{reduction = TRUE}.
#'   \item \code{Run.Dissimilarity()} — computes dissimilarities across methods.
#' }
#'
#' This provides a simple entry point for end-to-end setup and dissimilarity computation 
#' from raw single-cell data with minimal manual steps.
#'
#' @examples
#' \dontrun{
#' sc_res <- wrapper_scTypeEval(
#'   count_matrix = sc_counts,
#'   metadata = sc_metadata,
#'   ident = "celltype",
#'   sample = "sample_id",
#'   normalization.method = "Log1p",
#'   dissimilarity.method = c("WasserStein", "Pseudobulk:Euclidean"),
#'   reduction = TRUE,
#'   ndim = 30,
#'   verbose = TRUE
#' )
#'
#' # Access computed dissimilarities
#' names(sc_res@dissimilarity)
#' }
#'
#' @seealso 
#' \code{\link{create.scTypeEval}}, 
#' \code{\link{Run.ProcessingData}}, 
#' \code{\link{Run.Dissimilarity}}, 
#' \code{\link{Run.PCA}}, 
#' \code{\link{Run.HVG}}
#'
#' @export wrapper_scTypeEval

wrapper_scTypeEval <- function(scTypeEval = NULL,
                               count_matrix,
                               metadata,
                               ident,
                               sample,
                               aggregation = c("single-cell", "pseudobulk"),
                               gene.list = NULL,
                               reduction = TRUE,
                               ndim = 30,
                               black.list = NULL,
                               normalization.method = "Log1p",
                               dissimilarity.method = c("WasserStein", "Pseudobulk:Euclidean",
                                                        "Pseudobulk:Cosine", "Pseudobulk:Pearson",
                                                        "RecipClassif:Match", "RecipClassif:Score"),
                               min.samples = 5,
                               min.cells = 10,
                               ncores = 1,
                               bparam = NULL,
                               progressbar = FALSE,
                               verbose = TRUE){
   
   if(is.null(scTypeEval)){
      sc <- create.scTypeEval(matrix = count_matrix,
                              metadata = metadata,
                              active.ident = ident,
                              black.list = black.list)
   } else {
      sc <- scTypeEval
   }
   
   if(is.null(aggregation)){
      aggregation <- dissimilarity_methods[names(dissimilarity_methods) %in% dissimilarity.method] |>
         unique()
   }
   
   sc <- Run.ProcessingData(sc, 
                            ident = ident,
                            sample = sample,
                            aggregation = aggregation,
                            normalization.method = normalization.method,
                            min.samples = min.samples,
                            min.cells = min.cells,
                            verbose = verbose)
   
   if(is.null(gene.list)){
      if(verbose){message("Obtaining HVG from ", aggregation[1], "slot.")}
      sc <- Run.HVG(sc,
                    aggregation = aggregation[1],
                    ncores = ncores,
                    verbose = verbose)
   } else {
      sc <- add.GeneList(sc,
                         gene.list)
   }
   
   if(reduction){
      sc <- Run.PCA(sc,
                    ndim = ndim,
                    verbose = verbose)
   }
   
   for(m in dissimilarity.method){
      if(verbose){message(">.  Running ", m, "\n")}
      sc <- Run.Dissimilarity(sc,
                              reduction = reduction,
                              method = m,
                              ncores = ncores,
                              verbose = verbose)
   }
   
   return(sc)
   
}



#' @title Optimal Hierarchical Clustering Based on Consistency Metrics
#'
#' @description
#' An automated clustering procedure that iteratively splits cell populations
#' into subclusters and selects an optimal clustering resolution based on
#' cluster consistency metrics.
#'
#' This function performs a top-down hierarchical clustering strategy,
#' evaluates each split using one or more consistency metrics,
#' and retains only those splits that improve or maintain consistency above
#' user-defined thresholds. The final output is an \code{scTypeEval} object
#' annotated with an optimal clustering assignment.
#'
#' @param X Optional numeric matrix of features (cells as rows, features as columns). 
#'    Usually a low dimensional embeddings for all cells of datasets are expected.
#'   If \code{NULL}, the function preprocesses the provided \code{scTypeEval} object
#'   using \code{process_clustering()} to obtain a PCA embedding or filtered matrix.
#' @param scTypeEval A \code{scTypeEval} object containing raw or processed
#'   single-cell data and metadata.
#' @param sample Character string specifying the metadata column containing
#'   sample identifiers.
#' @param reduction Logical; if \code{TRUE}, PCA is performed during preprocessing
#'   when \code{X} is not provided (default: \code{TRUE}).
#' @param ndim Integer; number of principal components to retain when
#'   \code{reduction = TRUE} (default: \code{30}).
#' @param gene.list Optional named list of gene sets used for consistency
#'   computation. If \code{NULL}, highly variable genes (HVGs) are computed.
#' @param min.cells Integer; minimum number of cells required per cluster
#'   to be considered during consistency computation (default: \code{10}).
#' @param min.samples Integer; minimum number of samples required per cluster
#'   for consistency evaluation (default: \code{5}).
#' @param clustering_method Character string specifying the clustering algorithm
#'   used for splitting. Currently supported: \code{"kmeans"}.
#'   \code{"leiden"} is reserved for future implementation.
#' @param consistency_method Character vector specifying consistency metrics to
#'   evaluate cluster splits. Each entry must follow the format
#'   \code{"<metric> | <dissimilarity.method>"}.
#' @param hvg.ngenes Integer; number of highly variable genes to select when
#'   \code{gene.list = NULL} (default: \code{2000}).
#' @param normalization.method Character string specifying the normalization
#'   method used during preprocessing (default: \code{"Log1p"}).
#' @param ncores Integer; number of CPU cores for parallel execution
#'   (default: \code{1}).
#' @param nchild Integer; number of childs per cluster in each iteration when using kmeans (default: \code{2}).
#' @param max_Nclusters Integer; maximum number of clusters allowed before
#'   stopping the hierarchical splitting process (default: \code{100}).
#' @param max_iter Integer; maximum number of hierarchical splitting iterations
#'   (default: \code{100}).
#' @param nstart Integer; number of random starts for k-means clustering
#'   (default: \code{30}).
#' @param epsilon Numeric; tolerance allowing child cluster consistency to be
#'   slightly lower than the parent cluster consistency (default: \code{0.2}).
#' @param min.consistency Numeric; minimum consistency score required for a
#'   cluster to be retained in the final solution (default: \code{0.5}).
#' @param weight_by Character string specifying how to aggregate child
#'   consistency scores when evaluating a split:
#'   \itemize{
#'     \item \code{"none"} — unweighted mean
#'     \item \code{"cells"} — weighted by number of cells
#'     \item \code{"samples"} — weighted by number of samples
#'   }
#' @param verbose Logical; if \code{TRUE}, prints progress messages
#'   (default: \code{TRUE}).
#'
#' @return
#' An updated \code{scTypeEval} object containing:
#' \itemize{
#'   \item A new metadata column \code{optimal} with the selected clustering labels
#'   \item Intermediate hierarchical clustering assignments
#'   \item Consistency results stored in \code{scTypeEval@consistency}
#' }
#'
#' @details
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item Preprocesses the data if \code{X} is not provided
#'   \item Initializes all cells into a root cluster
#'   \item Iteratively splits clusters into subclusters
#'   \item Computes consistency metrics for each split using
#'         \code{compute_consistency()}
#'   \item Retains splits that satisfy consistency improvement and threshold
#'         criteria
#'   \item Reverts clusters failing \code{min.consistency} to their parent cluster
#' }
#'
#' This strategy identifies a data-driven clustering resolution without
#' requiring the number of clusters to be specified a priori.
#'
#' @examples
#' \dontrun{
#' scTypeEval <- optimal_clustering(
#'   scTypeEval = sc_obj,
#'   ident = "celltype",
#'   sample = "sample_id",
#'   consistency_method = c(
#'     "silhouette | RecipClassif:Match",
#'     "2label.silhouette | Pseudobulk:Cosine"
#'   ),
#'   min.consistency = 0.6,
#'   verbose = TRUE
#' )
#'
#' table(scTypeEval@metadata$optimal)
#' }
#'
#' @seealso
#' \code{\link{process_clustering}},
#' \code{\link{compute_consistency}},
#' \code{\link{get.clusters}}
#'
#' @export get.optimal_clustering


get.optimal_clustering <- function(X = NULL,
                                   scTypeEval,
                                   sample = "sample",
                                   reduction = TRUE,
                                   ndim = 30,
                                   gene.list = NULL,
                                   min.cells = 10,
                                   min.samples = 5,
                                   clustering_method = c("kmeans", "louvain", "leiden"),
                                   consistency_method = c("silhouette | RecipClassif:Match",
                                                          "2label.silhouette | Pseudobulk:Cosine"),
                                   hvg.ngenes = 2000,
                                   normalization.method = "Log1p",
                                   ncores = 1,
                                   nchild = 2,
                                   max_Nclusters = 100,
                                   max_iter = 100,
                                   nstart = 30,
                                   epsilon = 0.2,
                                   min.consistency = 0.5,
                                   weight_by = c("none", "cells", "samples" ),
                                   verbose = TRUE) {
   
   clustering_method <- clustering_method[1] |> tolower()
   # install igraph if not done
   if(clustering_method %in% c("louvain", "leiden") && !requireNamespace("igraph", quietly = TRUE)){
      message("igraph package needed for ", clustering_method., ". Installing igraph...\n")
      install.packages("igraph")
   }
   weight_by <- weight_by[1] |> tolower()
   
   # 1. preprocess data if X not provided
   if(is.null(X)){
      if(verbose) message("X not provided, processing scTypeEval object")
      ret <- process_clustering(scTypeEval,
                                sample = sample,
                                reduction = reduction,
                                ndim = ndim,
                                gene.list = gene.list,
                                hvg.ngenes = hvg.ngenes,
                                normalization.method = normalization.method,
                                ncores = ncores,
                                verbose = verbose)
      X <- t(ret$X)
      if(is.null(gene.list)) gene.list <- ret$gene.list
   }
   
   scTypeEval@metadata$.tmp <- "root"
   allcells <- rownames(X)
   ct_follow <- "root"
   cons.list <- list()
   parent_score <- c("root" = 0)
   iter <- 0
   
   while(dplyr::n_distinct(scTypeEval@metadata$.tmp) < max_Nclusters && iter < max_iter){
      
      if(length(ct_follow) == 0) {break}
      child_vector <- c()
      to_follow <- c()
      iter <- iter + 1
      clus_col <- paste(clustering_method, iter, sep = "_")
      for(cu in ct_follow){
         
         cells <- allcells[scTypeEval@metadata[[".tmp"]] == cu]
         
         resolution0 <- pick_resolution(length(cells), nchild)
         resolution <- resolution0
         resolution_iter <- 0 # max iteration trials for louvan/leiden resolution
         nclus <- 1
         while(nclus < 2 && resolution_iter <= 10){
            set.seed(22)
            cl <- get.clusters(X[cells, , drop=FALSE],
                               clustering_method = clustering_method,
                               nclusters = nchild,
                               nstart = nstart,
                               resolution = resolution,
                               ncores = ncores)
            
            nclus <- dplyr::n_distinct(cl)
            resolution_iter <- resolution_iter + 1
            resolution <- resolution0 + resolution0*resolution_iter
         }
         
         if(nclus<2) {next}
         
         scTypeEval@metadata[cells, ".tmp"] <- paste(scTypeEval@metadata[cells, ".tmp"], cl, sep = ".")
         if(verbose){cat("Computing consistency for", clus_col)}
         
         suppressMessages(
            {
               cons <- compute_consistency(scTypeEval,
                                           ident = ".tmp",
                                           sample = sample,
                                           gene.list = gene.list,
                                           consistency_method = consistency_method,
                                           min.samples = min.samples,
                                           min.cells = min.cells,
                                           ncores = ncores,
                                           verbose = F)
            })
         
         # aggregate child consistency
         child_celltypes <- unique(scTypeEval@metadata[cells, ".tmp"])
         cons <- cons |>
            dplyr::filter(celltype %in% child_celltypes) 
         child_cons <- cons |>
            dplyr::filter(celltype %in% child_celltypes) |>
            dplyr::pull(product, name = celltype)
         if(weight_by == "cells"){
            sizes <- table(scTypeEval@metadata[cells, ".tmp"])
            child_score <- weighted.mean(child_cons[child_celltypes], sizes[child_celltypes])
         } else if(weight_by == "samples"){
            nsamp <- scTypeEval@metadata[cells,] |>
               dplyr::group_by(.data[[".tmp"]]) |>
               dplyr::summarise(n = dplyr::n_distinct(.data[[sample]]))
            weights <- nsamp$n[match(names(child_cons), nsamp$.tmp)]
            child_score <- weighted.mean(child_cons, weights)
         } else {
            child_score <- mean(child_cons, na.rm = TRUE)
         }
         
         cons.list[[paste(clus_col, cu, sep = "-")]] <- cons |>
            dplyr::mutate(ident = clus_col,
                          parent = cu)
         
         if(child_score >= (parent_score[cu] - epsilon) && child_score >= min.consistency){
            to_follow <- c(to_follow, child_celltypes)
         }
         child_score_pass <- cons |>
            dplyr::filter(celltype %in% to_follow) |>
            dplyr::pull(product, name = celltype)
         child_vector <- c(child_vector, child_score_pass)
      }
      ct_follow <- to_follow
      parent_score <- child_vector
      scTypeEval@metadata[[clus_col]] <-  scTypeEval@metadata[[".tmp"]]
   }
   
   # get optimal clustering
   allcons <- do.call(rbind, cons.list)
   
   nfail <- 1
   iter_fail <- 0
   scTypeEval@metadata$optimal0 <- scTypeEval@metadata$.tmp
   
   while(length(nfail) > 0 && iter_fail < 500){
      iter_fail <- iter_fail + 1 # safety break
      # get clusters do not passing the threshold of min.consistency
      # revert all childs to parent, even if only one child is poor consistent
      nfail_parent <- allcons |>
         dplyr::filter(product < min.consistency ) |>
         dplyr::filter(celltype %in% unique(scTypeEval@metadata$optimal0)) |>
         dplyr::pull(parent)
      nfail <- allcons |>
         dplyr::filter(parent %in% nfail_parent) |>
         dplyr::pull(celltype)
      
      if(length(nfail) == 0) break  # stop if none fail
      
      # for clustering under the threshold, use parent clustering
      scTypeEval@metadata <- 
         scTypeEval@metadata |>
         dplyr::rowwise() |>
         dplyr::mutate(optimal0 = ifelse(optimal0 %in% nfail,
                                         allcons |> dplyr::filter(celltype == optimal0) |> dplyr::pull(parent),
                                         optimal0)
         ) |>
         dplyr::ungroup()
   }
   
   scTypeEval@metadata <- 
      scTypeEval@metadata |>
      dplyr::mutate(optimal = factor(optimal0),
                    optimal = paste0("C", as.integer(optimal)))
   
   # add consistency output
   scTypeEval@consistency <- c(scTypeEval@consistency, cons.list)
   
   return(scTypeEval)
}
