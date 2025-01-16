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


split.matrix <- function(mat,
                         ident,
                         sample = NULL,
                         min.cells = 10,
                         bparam = BiocParallel::SerialParam()){
   if (length(sample) != ncol(mat)) {
      stop("The length of 'sample' must match the number of columns in 'mat'")
   }
   
   # Get the unique sample groups
   unique_samples <- levels(sample)
   
   # Initialize an empty list to store matrices
   split_list <- list()
   
   # Loop through each unique sample and subset the matrix
   split_list <- BiocParallel::bplapply(unique_samples,
                                        BPPARAM = bparam,
                                        function(s){
                                           new.mat <- mat[, sample == s, drop = FALSE]
                                           new.ident <- ident[sample == s]
                                           
                                           ret <- valid.sc(mat = new.mat,
                                                           ident = new.ident,
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
                              bparam = BiocParallel::SerialParam(),
                              min.samples = 5,
                              min.cells = 10,
                              KNNGraph_k = 5){
   
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
                                        
                                        
                                        # compute internal validation metrics
                                        con <- calculate_IntVal_metric(mat = red.mat@matrix,
                                                                       norm.mat = norm.mat,
                                                                       metrics = IntVal.metric,
                                                                       distance.method = distance.method,
                                                                       ident = red.mat@ident,
                                                                       KNNGraph_k = KNNGraph_k)
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
   if (inverse) {
      # For values where higher means worse clustering, scale to [0, 1] (inverted)
      return(1 - ((value - min_value) / (max_value - min_value)))
   } else {
      # For values where higher means better clustering, scale to [0, 1]
      return((value - min_value) / (max_value - min_value))
   }
}

normalize_metric <- function(value, metric) {
   # Normalize all metrics to the range [0, 1]
   scaled_metric <- switch(metric,
                           
                           # remove negative values
                           "silhouette" = rm0(value),
                           # remove negative values
                           "modularity" = rm0(value),
                           "modularity_pct" = minmax_norm(value,
                                                      min_value = -1,
                                                      max_value = 1,
                                                      inverse = FALSE),
                           # Ward, neighborhood purity and BestHit are already in [0, 1],
                           # no need to scale
                           "NeighborhoodPurity" = rm0(value),
                           "ward.PropMatch" = rm0(value),
                           "Leiden.PropMatch" = rm0(value),
                           "BestHit-Mutual.Score" = value,
                           "BestHit-Mutual.Match" = rm0(value),
                           
                           # Normalize inertia to [0, 1] (higher is worse, hence inverse = TRUE)
                           "inertia" = minmax_norm(value,
                                                   min_value = min(value),
                                                   max_value = max(value),
                                                   inverse = TRUE),
                           
                           # Normalize Xie-Beni to [0, 1] (higher is worse, hence inverse = TRUE)
                           "Xie-Beni" = minmax_norm(value,
                                                    min_value = min(value),
                                                    max_value = max(value),
                                                    inverse = TRUE),
                           
                           # Normalize S_Dbw to [0, 1] (higher is worse, hence inverse = TRUE)
                           "S_Dbw" = minmax_norm(value,
                                                 min_value = min(value),
                                                 max_value = max(value),
                                                 inverse = TRUE),
                           
                           # Normalize I to [0, 1] (higher is better, inverse = FALSE)
                           "I" = minmax_norm(value,
                                             min_value = min(value),
                                             max_value = max(value),
                                             inverse = FALSE)
   )
   
   return(scaled_metric)
}



get.PCA <- function(mat,
                    ident,
                    normalization.method,
                    data.type,
                    sample = NULL,
                    min.samples = 5,
                    min.cells = 10,
                    bparam = BiocParallel::SerialParam(),
                    ndim = 50)
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
                                        
                                        # if ncol or nrow is below given ndim use this number
                                        ndim <- min(dim(norm.mat), ndim)
                                        
                                        # compute PCA
                                        pr <- stats::prcomp(Matrix::t(norm.mat),
                                                            rank. = ndim)
                                        
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
