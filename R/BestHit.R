mutual_method <-  c("Mutual.Score", "Mutual.Match")

# code to perform a best-hit classifier based on singleR
singleR.helper <- function(test,
                           ref,
                           bparam = BiocParallel::SerialParam()){
   pred <- SingleR::SingleR(
      test = test,
      ref = ref,
      labels = ref$ident,
      de.method = "classic",
      BPPARAM = bparam
   ) |> 
      as.data.frame() 
   return(pred)
}

singleR.score.tidy <- function(pred,
                               .true){
   pred <- pred |>
      dplyr::select(dplyr::starts_with("score")) |>  
      dplyr::mutate(true = .true) |>
      tidyr::pivot_longer(-c(true),
                          names_to = "celltype",
                          values_to = "score") |>
      dplyr::filter(true == gsub("scores[.]", "", celltype)) |>
      dplyr::group_by(true) |>
      dplyr::summarize(score = mean(score))
   
   return(pred)
}

singleR.score <- function(pred1, true1,
                          pred2, true2){
   
   pred1 <- singleR.score.tidy(pred1, true1)
   pred2 <- singleR.score.tidy(pred2, true2)
   
   pred <- dplyr::inner_join(pred1, pred2, by = "true") |>
      dplyr::mutate(score = score.x * score.y) |>
      dplyr::select(true, score) |>
      dplyr::rename("celltype" = "true")
   
   return(pred)
}

singleR.match.tidy <- function(pred, .true){
   pred <- pred |>
      dplyr::select("pruned.labels") |>
      dplyr::mutate(true = .true,
                    score = ifelse(true == pruned.labels,
                                   1, 0))
   
   return(pred)
}


singleR.match <- function(pred1, pred2,
                          true1, true2){
   
   pred1 <- singleR.match.tidy(pred1, true1)
   pred2 <- singleR.match.tidy(pred2, true2)
   
   pred <- dplyr::inner_join(pred1, pred2, by = "true") |>
      dplyr::group_by(true) |>
      dplyr::mutate(score = score.x * score.y,
                    score = sum(score[score == 1])/dplyr::n()) |>
      dplyr::select(true, score) |>
      dplyr::rename("celltype" = "true")

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
   sce <- sce[,Matrix::colSums(SingleCellExperiment::counts(sce)) > 0] # Remove libraries with no counts.
   sce <- scuttle::logNormCounts(sce)
   return(sce)
}

bestHit.SingleR <- function(mat,
                            ident,
                            ident_GroundTruth = NULL,
                            sample,
                            method = "Mutual.Score",
                            data.type = "sc",
                            min.cells = 10,
                            min.samples = 5,
                            sep = "_",
                            ncores = 1,
                            bparam = NULL,
                            progressbar = TRUE){
   
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
      nc <- lapply(sce.listGT[c(a,b)],
                   function(x){length(unique(x$ident))})
      if(any(nc < 2)){ next }
      
      pred1 <- singleR.helper(sce.list[[a]],sce.listGT[[b]],bparam = param)
      pred2 <- singleR.helper(sce.list[[b]],sce.listGT[[a]], bparam = param)
      
      # original annotations to compare with predictions
      true1 <- sce.list[[a]]$ident
      true2 <- sce.list[[b]]$ident
      
      results <- BiocParallel::bplapply(method,
                                        BPPARAM = param,
                                        function(meth){
      
       r <- switch(meth,
                   "Mutual.Score" = singleR.score(pred1 = pred1,
                                                  pred2 = pred2,
                                                  true1 = true1,
                                                  true2 = true2),
                   "Mutual.Match" = singleR.match(pred1 = pred1,
                                                  pred2 = pred2,
                                                  true1 = true1,
                                                  true2 = true2),
                   stop(meth, " is not a supported Mutual Besthit method."))
       return(r)
                                        })
      names(results) <- method

      
      df.tmp[[paste(a, b, sep = "_vs_")]] <- results
   }
   
   # join results by method
   resli <- lapply(method, function(m)
   {
      res <- do.call(rbind, lapply(df.tmp, function(x) x[[m]]))
      # summarize results per cell type
      res <- res |>
         dplyr::group_by(celltype) |>
         dplyr::summarize(score = mean(score, na.rm = T))
      
      if(m == "Mutual.Match"){
         # normalize for BestHiT match
         nident <- sapply(mat.split, function(mat) {length(mat@ident)})
         # Generate all pairwise combinations of the samples
         sample_combinations <- combn(nident, 2)
         # Calculate pairwise probabilities (1 / (cell types in sample i * cell types in sample j))
         probabilities <- apply(sample_combinations, 2, function(pair) {
            1 / (pair[1] * pair[2])})
         # Average all pairwise probabilities together
         min.val <- mean(probabilities) 
         res <- res |>
            dplyr::mutate(score = minmax_norm(score,
                                              min_value = min.val,
                                              max_value = 1))
      }
      # extract vector with output
      res <-  res |>
         dplyr::pull(score, name = celltype)
      
      return(res)
   }
   )
   names(resli) <- method
   
   return(resli)
   
}
