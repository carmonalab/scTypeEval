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
                              min.cells = 10){
   
   mat <- get.matrix(mat,
                     data.type = data.type,
                     ident = ident,
                     sample = sample,
                     min.samples = min.samples,
                     min.cells = min.cells,
                     bparam = bparam)

   if(!is.list(mat)){
      mat <- list(mat)
   }
   
   consist <- BiocParallel::bplapply(mat,
                                     BPPARAM = bparam,
                                     function(red.mat){
                                        
                                        # normalized subetted matrix
                                        norm.mat <- Normalize_data(red.mat@matrix,
                                                                   method = normalization.method)
                                        
                                        
                                        # compute internal validation metrics
                                        con <- calculate_IntVal_metric(mat = red.mat@matrix,
                                                                       norm.mat = norm.mat,
                                                                       metrics = IntVal.metric,
                                                                       dist.method = distance.method,
                                                                       ident = red.mat@ident)
                                        return(con)
                                        
                                     })
   
   names(consist) <- names(mat)
   
   # combine consistencies if 1vsAll
   if(data.type == "pseudobulk_1vsall"){
      join <- list()
      for(cm in IntVal.metric){
         vec <- lapply(consist, function(m){
            m[[cm]][names(m[[cm]]) != "psblk"]
         }) %>% unlist()
         
         join[[cm]] <- vec
      }
      consist <- join
   }

   
   return(consist)
   
}


