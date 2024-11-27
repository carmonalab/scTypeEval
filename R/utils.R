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


consistency.helper <- function(mat,
                              ident,
                              normalization.method,
                              IntVal.metric,
                              data.type,
                              sample = NULL){
   
   mat <- get.matrix(mat,
                     data.type = data.type,
                     ident = ident,
                     sample = sample,
                     min.samples = min.samples,
                     min.cells = min.cells)
   
   
   
   if(!is.list(mat)){
      mat <- as.list(mat)
   }
   
   consist <- lapply(mat,
                     function(red.mat){
                        
                        # normalized subetted matrix
                        norm.mat <- Normalize_data(red.mat[["matrix"]],
                                                   method = normalization.method)
                        
                        # compute internal validation metrics
                        con <- calculate_IntVal_metric(mat = red.mat[["matrix"]],
                                                       norm.mat = norm.mat,
                                                       metrics = IntVal.metric,
                                                       dist.method = distance.method,
                                                       ident = red.mat[["ident"]])
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
   
   names(consist) <- IntVal.metric
   
   return(consist)
   
}


