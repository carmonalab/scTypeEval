# code to perform a best-hit classifier based on singleR

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
   
   res <- data.frame(matrix(nrow= 0, ncol = length(sce.list)+1))
   names(res) <- c("original", paste("sample", 1:(length(sce.list)-1), sep = "_"), "hit_proportion")
   for(a in names(sce.list)){
      df.tmp <- data.frame(original = sce.list[[a]]$ident)
      s <- 1
      for(b in names(sce.list)[names(sce.list) != a]){
         pred <- SingleR::SingleR(
            test = sce.list[[a]],
            ref = sce.list[[b]],
            labels = sce.list[[b]]$ident,
            BPPARAM = param
         ) |> 
            as.data.frame() |>
            dplyr::mutate(pruned.labels = ifelse(is.na(pruned.labels),
                                                       "na",
                                                       pruned.labels)
                          ) |>
           dplyr::pull(pruned.labels)
         
         df.tmp[[paste("sample", s, sep = "_")]] <- pred
         s <- s + 1
      }
      df.tmp <- df.tmp |>
         dplyr::rowwise() |>
         dplyr::mutate(
            hit_proportion = sum(c_across(starts_with("sample")) == original, na.rm = TRUE)/(length(sce.list)-1)
         ) |>
         ungroup()
      res <- rbind(res, df.tmp)
   }
   
   # summarize results per cell type
   res <- res |>
      dplyr::group_by(original) |>
      dplyr::summarize(hit_proportion = mean(hit_proportion)) |>
      dplyr::pull(hit_proportion, name = original)
   
   return(res)
   
}
