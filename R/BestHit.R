# code to perform a best-hit classifier based on singleR
singleR.helper <- function(test,
                           ref,
                           bparam = BiocParallel::SerialParam()){
   pred <- SingleR::SingleR(
      test = test,
      ref = ref,
      labels = ref$ident,
      BPPARAM = bparam
   ) |> 
      as.data.frame() |> 
      dplyr::select(dplyr::starts_with("score"))
   pred$cellid <- rownames(pred)
   pred <- pred |>  
      dplyr::mutate(true = test$ident) |>
      tidyr::pivot_longer(-c(cellid, true),
                          names_to = "celltype",
                          values_to = "score") |>
      dplyr::filter(true == gsub("scores[.]", "", celltype)) |>
      dplyr::group_by(true) |>
      dplyr::summarize(score = mean(score))
   
   return(pred)
   
}



bestHit.SingleR <- function(mat,
                            ident,
                            sample,
                            min.cells = 10,
                            ncores = 1,
                            bparam = NULL,
                            progressbar = TRUE){
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   mat.split <- split.matrix(mat,
                             ident = ident,
                             sample = sample,
                             min.cells = min.cells,
                             bparam = param)
   
   sce.list <- BiocParallel::bplapply(mat.split,
                                      BPPARAM = param,
                                      function(m){
                                         sce <- SingleCellExperiment::SingleCellExperiment(
                                            assays = list(counts = m@matrix),
                                            colData = data.frame(ident = m@ident)
                                         )
                                         sce <- sce[,colSums(SingleCellExperiment::counts(sce)) > 0] # Remove libraries with no counts.
                                         sce <- scuttle::logNormCounts(sce)
                                         return(sce)
                                      })
   names_list <- names(sce.list)
   # Create combinations and ensure columns are not factors
   combis <- expand.grid(names_list, names_list, stringsAsFactors = FALSE) |>
      dplyr::filter(Var1 < Var2)
   
   df.tmp <- list()
   
   for(i in 1:nrow(combis)){
      a <- combis[i,1]
      b <- combis[i,2]
      
      # if reference consist of only one cell, singleR do not return score
      # but all cells are classified as this unique cell with score of 1
      nc <- lapply(sce.list[c(a,b)],
                   function(x){length(unique(x$ident))})
      if(any(nc < 2)){ next }
      
      pred1 <- singleR.helper(sce.list[[a]],sce.list[[b]],bparam = param)
      pred2 <- singleR.helper(sce.list[[b]],sce.list[[a]], bparam = param)
      
      pred <- dplyr::inner_join(pred1, pred2, by = "true") |>
         dplyr::mutate(score = score.x * score.y,
                       comp = paste(a, b, sep = "_vs_")) |>
         dplyr::select(true, score, comp) |>
         dplyr::rename("celltype" = "true")
      
      df.tmp[[paste(a, b, sep = "_vs_")]] <- pred
   }
   
   res <- do.call(rbind, df.tmp)
   # summarize results per cell type
   res <- res |>
      dplyr::group_by(celltype) |>
      dplyr::summarize(score = mean(score)) |>
      dplyr::pull(score, name = celltype)
   
   return(res)
   
}
