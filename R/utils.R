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
   
   # adjust for GloScope we require still sc, but filtered
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
   
   # filter for mats with at least 2 cell types
   mats <- Filter(function(m) length(levels(m@ident)) > 1, mats)
   
   consist <- BiocParallel::bplapply(mats,
                                     BPPARAM = bparam,
                                     function(red.mat){
                                        # normalized subetted matrix
                                        norm.mat <- Normalize_data(red.mat@matrix,
                                                                   method = normalization.method)
                                        dist <- NULL
                                        nident <- red.mat@ident
                                        if(pca){
                                           if(verbose){message("Computing PCA space\n")}
                                           pr <- custom_prcomp(norm.mat, ndim)
                                           mat <- NULL
                                           
                                           if(data.type == "GloScope"){
                                              if(verbose){message("Computing GloScope divergencies\n")}
                                              suppressWarnings(
                                                 {
                                                    dist <- GloScope::gloscope(embedding_matrix = pr$x,
                                                                               cell_sample_ids = red.mat@ident,
                                                                               dens = "KNN",
                                                                               dist_mat = distance.method,
                                                                               k = min.cells-1
                                                    )
                                                    nident <- sapply(colnames(dist),
                                                                     function(x){
                                                                        strsplit(x, "_")[[1]][2]
                                                                     }) |>
                                                       factor()
                                                 })
                                              norm.mat <- dist
                                              dist <- as.dist(dist)
                                           } else {
                                              norm.mat <- t(pr$x)
                                           }
                                        } else{
                                           mat <- red.mat@matrix
                                        }
                                        
                                        
                                        # compute internal validation metrics
                                        con <- calculate_IntVal_metric(mat = mat,
                                                                       norm.mat = norm.mat,
                                                                       dist = dist,
                                                                       metrics = IntVal.metric,
                                                                       distance.method = distance.method,
                                                                       ident = nident,
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
                    normalization.method = "Log1p",
                    aggregation = c("single-cell"),
                    sample = NULL,
                    gene.list,
                    min.samples = 5,
                    min.cells = 10,
                    black.list = NULL,
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
                          ndim, 
                          verbose = TRUE){
   # if ncol or nrow is below given ndim use this number
   ndim <- min(dim(norm.mat)-1, ndim)
   
   if(verbose){message("Returning ", ndim, "dimensions for PCA")}
   
   # compute PCA
   pr <- irlba::prcomp_irlba(Matrix::t(norm.mat),
                             n = ndim,
                             center = TRUE,
                             scale. = FALSE)
   return(pr)
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
                                         layer = "counts")
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



sample_variable_length_combinations <- function(elements,
                                                min_k = 2,
                                                max_k = NULL,
                                                num_samples = 30,
                                                max_iter = 1000,
                                                seed = 22) {
   if (is.null(seed)) {
      stop("Seed value missing")
   }
   set.seed(seed)
   
   if(is.null(max_k)){
      max_k <- length(elements) - 1
   }
   
   unique_combinations <- list()
   count <- 0
   iter <- 0
   
   # A simple way to sample k values (could be weighted if desired)
   possible_k_values <- min_k:max_k
   
   prop0 <- rep(1, length(elements))
   names(prop0) <- elements
   
   while (count < num_samples) {
      # 1. Randomly pick a size 'k' for the combination
      current_k <- sample(possible_k_values, 1)
      
      # Ensure k is not larger than elements available
      if (current_k > length(elements)) next
      
      
      # 2. Sample 'current_k' elements
      sampled_elements <- sample(elements,
                                 current_k,
                                 replace = FALSE)
      
      # 3. Sort to make it a canonical representation of a combination
      sorted_combination <- sort(sampled_elements)
      current_combo_string <- paste(sorted_combination, collapse = "-")
      
      if (!current_combo_string %in% names(unique_combinations)) {
         unique_combinations[[current_combo_string]] <- sorted_combination
         count <- count + 1
      }
      
      iter <- iter + 1
      
      if(iter >= max_iter){
         break
      }
   }
   
   return(unique_combinations)
}


get.MDS <- function(mat,
                    ident,
                    normalization.method = "Log1p",
                    distance.method = "euclidean",
                    data.type = "sc",
                    sample = NULL,
                    min.samples = 5,
                    min.cells = 10,
                    bparam = BiocParallel::SerialParam(),
                    ndim = 30,
                    pca = FALSE,
                    pca.dim = 30,
                    verbose = TRUE)
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
   
   mdss <- BiocParallel::bplapply(mats,
                                  BPPARAM = bparam,
                                  function(red.mat){
                                     # normalized subetted matrix
                                     norm.mat <- Normalize_data(red.mat@matrix,
                                                                method = normalization.method)
                                     nident <- red.mat@ident
                                     if(pca){
                                        if(verbose){message("Computing PCA space\n")}
                                        pr <- custom_prcomp(norm.mat, pca.dim)
                                        mat <- NULL
                                        
                                        if(data.type == "GloScope"){
                                           if(verbose){message("Computing GloScope divergencies\n")}
                                           suppressWarnings(
                                              {
                                                 dist <- GloScope::gloscope(embedding_matrix = pr$x,
                                                                            cell_sample_ids = red.mat@ident,
                                                                            dens = "KNN",
                                                                            dist_mat = distance.method,
                                                                            k = min.cells-1
                                                 )
                                                 nident <- sapply(colnames(dist),
                                                                  function(x){
                                                                     strsplit(x, "_")[[1]][2]
                                                                  }) |>
                                                    factor()
                                              })
                                           norm.mat <- dist
                                           dist <- as.dist(dist)
                                        } else {
                                           norm.mat <- t(pr$x)
                                        }
                                     } else{
                                        mat <- red.mat@matrix
                                     }
                                     # get distances
                                     if(data.type != "GloScope"){
                                        dist <- get.distance(mat = mat,
                                                             norm.mat = norm.mat,
                                                             distance.method = distance.method)
                                     }
                                     
                                     ndim <- min(attr(dist, "Size") - 1, ndim)
                                     
                                     # comute MDS scale
                                     mds <- stats::cmdscale(dist,
                                                            k = ndim,
                                                            eig = TRUE)
                                     
                                     rr <- methods::new("DimRed",
                                                        embeddings = mds$points,
                                                        feature.loadings = matrix(nrow = 0, ncol = 0),  # not applicable for MDS,
                                                        gene.list = "tmp",
                                                        black.list = "tmp",
                                                        data.type = as.character(data.type),
                                                        ident = nident,
                                                        key = "MDS")
                                     
                                     return(rr)
                                     
                                  })
   
   names(mdss) <- names(mats)
   
   return(mdss)
   
}


.general_filtering <- function(mat, # Mat_ident object
                              black.list = NULL,
                              gene.list = NULL){
   norm.mat <- mat@matrix
   # remove blacked listed genes
   if(!is.null(black.list) && verbose){message("   Filtering out black listed genes... \n")}
   norm.mat <- norm.mat[!rownames(norm.mat) %in% black.list,]
   
   # keep only gene list features
   if(verbose){message("   Filtering gene list... \n")}
   norm.mat <- norm.mat[rownames(norm.mat) %in% gene.list,]
   
   # remove rows or columns with only 0
   mat@matrix <- norm.mat
   if(verbose){message("   Filtering empty rows and cols... \n")}
   mat <- filter_empty(mat)
   return(mat)
}
