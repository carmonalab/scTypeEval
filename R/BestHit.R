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

get.MutualMatrix <- function(mat,
                             ident,
                             sample,
                             data.type = "sc",
                             min.cells = 10,
                             min.samples = 5,
                             sep = "_",
                             bparam = BiocParallel::SerialParam()){
   # Combined grouping
   groups <- factor(paste(sample, ident, sep = sep))
   
   valid_indices <- valid.indices(groups = groups,
                                  ident = ident,
                                  sample = sample,
                                  min.samples = min.samples,
                                  min.cells = min.cells,
                                  sep = sep)

   mat <- mat[, valid_indices, drop = FALSE]
   
   mats <- split.matrix(mat,
                        ident = ident[valid_indices],
                        sample = sample[valid_indices],
                        min.cells = min.cells,
                        bparam = bparam)
   
   mat.split <- lapply(mats,
                       function(mat){
                          get.matrix(mat@matrix,
                                     data.type = data.type,
                                     ident = mat@ident,
                                     sample = NULL, # now is pseudobulk per cell type per invidividual sample
                                     min.cells = min.cells,
                                     bparam = bparam)
                       })
   
   return(mat.split)
}

get.SCE <-  function(m){
   sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = m@matrix),
      colData = data.frame(ident = m@ident)
   )
   sce <- sce[,colSums(SingleCellExperiment::counts(sce)) > 0] # Remove libraries with no counts.
   sce <- scuttle::logNormCounts(sce)
   return(sce)
}

bestHit.SingleR <- function(mat,
                            ident,
                            ident_GroundTruth = NULL,
                            sample,
                            data.type = "sc",
                            min.cells = 10,
                            min.samples = 5,
                            sep = "_",
                            ncores = 1,
                            bparam = NULL,
                            progressbar = FALSE){
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   mat.split <- get.MutualMatrix(mat = mat,
                                 ident = ident,
                                 sample = sample,
                                 data.type = data.type,
                                 min.cells = min.cells,
                                 min.samples = min.samples,
                                 sep = sep,
                                 bparam = param)
   
   sce.list <- BiocParallel::bplapply(mat.split,
                                      BPPARAM = param,
                                      get.SCE)
   
   if(is.null(ident_GroundTruth)){
      sce.listGT <- sce.list
   } else {
      mat.split <- get.MutualMatrix(mat = mat,
                                    ident = ident_GroundTruth,
                                    sample = sample,
                                    data.type = data.type,
                                    min.cells = min.cells,
                                    min.samples = min.samples,
                                    sep = sep,
                                    bparam = param)
      
      sce.listGT <- BiocParallel::bplapply(mat.split,
                                           BPPARAM = param,
                                           get.SCE)
   }
   
   
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
      
      pred1 <- singleR.helper(sce.list[[a]],sce.listGT[[b]],bparam = param)
      pred2 <- singleR.helper(sce.list[[b]],sce.listGT[[a]], bparam = param)
      
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
