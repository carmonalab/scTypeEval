set_parallel_params <- function(ncores = NULL,
                                bparam = NULL,
                                progressbar = TRUE)
{
   if (is.null(ncores)) {
      ncores <- 1
   }
   
   if (ncores > parallel::detectCores()) {
      ncores <- parallel::detectCores()
      message("Using more cores available in this computer, reducing number of cores to ",
              ncores)
   }
   
   # set parallelization parameters
   if (is.null(bparam)) {
      if (ncores > 1) {
         if(.Platform$OS.type == "windows"){
            param <- BiocParallel::SerialParam(progressbar = progressbar)
            warning("Parallel processing not possible for Windows operating system.\n")
         } else {
            param <- BiocParallel::MulticoreParam(workers =  ncores,
                                                  progressbar = progressbar)
         }
      } else {
         param <- BiocParallel::SerialParam(progressbar = progressbar)
      }
   } else {
      param <- bparam
   }
   return(param)
}


split_matrix <- function(mat,
                         ident,
                         sample = NULL,
                         min.cells = 10,
                         bparam = BiocParallel::SerialParam()){
   if (length(sample) != ncol(mat)) {
      stop("The length of 'sample' must match the number of columns in 'mat'")
   }
   
   # Get the unique sample groups
   unique_samples <- unique(sample)
   
   # Initialize an empty list to store matrices
   split_list <- list()
   
   # Loop through each unique sample and subset the matrix
   split_list <- BiocParallel::bplapply(unique_samples,
                                        BPPARAM = bparam,
                                        function(s){
                                           new.mat <- mat[, sample == s, drop = FALSE]
                                           new.ident <- ident[sample == s]
                                           
                                           ret <- valid.sc(mat = new.mat,
                                                           ident = factor(new.ident),
                                                           min.cells = min.cells)
                                           
                                           return(ret)
                                        })
   names(split_list) <- as.character(unique_samples)
   
   return(split_list)
}


consistency.helper <- function(mat,
                              ident,
                              normalization.method,
                              distance.method = "euclidean",
                              IntVal.metric,
                              data.type,
                              sample = NULL,
                              pca = FALSE,
                              ndim = 30,
                              bparam = BiocParallel::SerialParam(),
                              min.samples = 5,
                              min.cells = 10,
                              KNNGraph_k = 5,
                              verbose = TRUE){
   
   mats <- get.matrix(mat,
                     data.type = data.type,
                     ident = ident,
                     sample = sample,
                     min.samples = min.samples,
                     min.cells = min.cells,
                     bparam = bparam)

   if(!is.list(mats)){
      mats <- list(mats)
   }
   
   consist <- BiocParallel::bplapply(mats,
                                     BPPARAM = bparam,
                                     function(red.mat){
                                        # normalized subetted matrix
                                        norm.mat <- Normalize_data(red.mat@matrix,
                                                                   method = normalization.method)
                                        if(!pca){
                                           mat <- red.mat@matrix
                                        } else {
                                           pr <- custom_prcomp(norm.mat, ndim)
                                           mat <- NULL
                                           norm.mat <- t(pr$x)
                                        }
                                        
                                        # compute internal validation metrics
                                        con <- calculate_IntVal_metric(mat = mat,
                                                                       norm.mat = norm.mat,
                                                                       metrics = IntVal.metric,
                                                                       distance.method = distance.method,
                                                                       ident = red.mat@ident,
                                                                       KNNGraph_k = KNNGraph_k,
                                                                       verbose = verbose)
                                        return(con)
                                        
                                     })
   
   names(consist) <- names(mats)
   
   # combine consistencies if 1vsAll
   if(data.type == "pseudobulk_1vsall"){
      join <- list()
      for(cm in IntVal.metric){
         vec <- lapply(consist, function(m){
            m[[cm]][names(m[[cm]]) != "psblk"]
         })
         vec <- sapply(vec, unname)
         join[[cm]] <- vec
      }
      consist <- list(join)
   } else{
      unlist(consist, recursive = F)
   }

   
   return(consist)
   
}

rm0 <- function(value){
   value[value < 0] <- 0
   return(value)
}


# Function to scale any value to the range [0, 1]
minmax_norm <- function(value, min_value, max_value, inverse = FALSE) {
   norm_value <- (value - min_value) / (max_value - min_value)
   norm_value <- pmax(0, pmin(1, norm_value))  # Clip values to [0,1]
   
   if (inverse) {
      return(1 - norm_value)
   } else {
      return(norm_value)
   }
}

normalize_metric <- function(value, metric) {
   # Define metric groups
   rm0_metrics <- c(
      "silhouette", "modularity", "modularity_pct",
      "NeighborhoodPurity", "ward.PropMatch", "Leiden.PropMatch",
      "Orbital.centroid", "Orbital.medoid",
      "ward.ARI", "Leiden.ARI",
      "GraphConnectivity", "ward.NMI", "Leiden.NMI",
      "BestHit-Mutual.Match"
   )
   
   minmax_inverse_metrics <- c("inertia", "Xie-Beni", "S.Dbw")
   minmax_direct_metrics <- c("I")
   
   # Normalize based on the metric group
   if (metric %in% rm0_metrics) {
      return(rm0(value))
   } else if (metric == "BestHit-Mutual.Score") {
      return(value)
   } else if (metric %in% minmax_inverse_metrics) {
      return(minmax_norm(value, min_value = min(value), max_value = max(value), inverse = TRUE))
   } else if (metric %in% minmax_direct_metrics) {
      return(minmax_norm(value, min_value = min(value), max_value = max(value), inverse = FALSE))
   } else {
      stop("Unknown metric: ", metric)
   }
}

get.PCA <- function(mat,
                    ident,
                    normalization.method,
                    data.type,
                    sample = NULL,
                    min.samples = 5,
                    min.cells = 10,
                    bparam = BiocParallel::SerialParam(),
                    ndim = 30)
{
   mats <- get.matrix(mat,
                      data.type = data.type,
                      ident = ident,
                      sample = sample,
                      min.samples = min.samples,
                      min.cells = min.cells,
                      bparam = bparam)
   
   if(!is.list(mats)){
      mats <- list(mats)
   }
   
   pcas <- BiocParallel::bplapply(mats,
                                     BPPARAM = bparam,
                                     function(red.mat){
                                        
                                        # normalized subetted matrix
                                        norm.mat <- Normalize_data(red.mat@matrix,
                                                                   method = normalization.method)
                                        
                                        pr <- custom_prcomp(norm.mat, ndim)
                                        
                                        rr <- methods::new("DimRed",
                                                           embeddings = pr$x,
                                                           feature.loadings = pr$rotation,
                                                           gene.list = "tmp",
                                                           black.list = "tmp",
                                                           data.type = as.character(data.type),
                                                           ident = red.mat@ident,
                                                           key = "PCA")
                                        
                                        return(rr)
                                        
                                     })
   
   names(pcas) <- names(mats)
   
   return(pcas)
   
}

custom_prcomp <- function(norm.mat,
                          ndim){
   # if ncol or nrow is below given ndim use this number
   ndim <- min(dim(norm.mat), ndim)
   
   # compute PCA
   pr <- stats::prcomp(Matrix::t(norm.mat),
                       rank. = ndim)
   return(pr)
}

# function to purge sample and annotation names
purge_label <- function(label){
   gsub(" |_|[+]", ".", label)
}


# function to compute variance captured in PC components

var_PCA <- function(pca_embeddings){
   # Compute variance per PC
   n_samples <- nrow(pca_embeddings)  # Number of samples
   variance_per_pc <- colSums(pca_embeddings^2) / (n_samples - 1)
   
   # Compute proportion of variance explained
   total_variance <- sum(variance_per_pc)
   variance_explained <- variance_per_pc / total_variance

   return(variance_explained)
}

# function to load single-cell objects
load_sc <- function(path,
                    split = TRUE) {
   
   if (!file.exists(path)) stop("File does not exist: ", path)
   
   # Initialize
   object <- NULL
   counts <- NULL
   metadata <- NULL
   
   if (grepl("rds$", path)) {
      object <- readRDS(path)
      
      if (inherits(object, "Seurat")) {
         rcounts <- Seurat::GetAssayData(object,
                                         assay = "RNA",
                                         slot = "counts")
         counts <- as(rcounts, "dgCMatrix")
         metadata <- as.data.frame(object@meta.data)
         
      } else if (inherits(object, "SingleCellExperiment")) {
         rcounts <- SummarizedExperiment::assay(object, "counts")
         counts <- as(rcounts, "dgCMatrix")
         metadata <- as.data.frame(SummarizedExperiment::colData(object))
         
      } else {
         stop("Unsupported .rds object type: ", class(object))
      }
      
   } else if (grepl("h5ad$", path)) {
      if (!requireNamespace("anndata", quietly = TRUE)) {
         stop("The 'anndata' package is required to read .h5ad files. Please install it.")
      }
      
      object <- anndata::read_h5ad(path)
      
      if (!"counts" %in% names(object$layers)) {
         stop("The .h5ad file does not contain a 'counts' layer.")
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

