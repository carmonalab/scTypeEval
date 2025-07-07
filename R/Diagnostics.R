hierarchy.helper <- function(mat,
                             ident,
                             normalization.method,
                             distance.method = "euclidean",
                             hierarchy.method = "ward.D2",
                             data.type = "pseudobulk",
                             sample = NULL,
                             pca = FALSE,
                             ndim = 30,
                             min.samples = 5,
                             min.cells = 10,
                             verbose = TRUE){
   
   mat <- get.matrix(mat,
                     data.type = data.type,
                     ident = ident,
                     sample = sample,
                     min.samples = min.samples,
                     min.cells = min.cells)
   
   # normalized subetted matrix
   norm.mat <- Normalize_data(red.mat@matrix,
                              method = normalization.method)
   nident <- red.mat@ident
   if(pca){
      if(verbose){message("Computing PCA space\n")}
      pr <- custom_prcomp(norm.mat, pca.dims)
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
         dist <- as.dist(dist)
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
   
   # change labels
   attr(dist, "Labels") <- nident
   
   hclust_result <- stats::hclust(dist,
                                  method = hierarchy.method)
   
   return(hclust_result)
   
}


nn.helper <- function(mat,
                      ident,
                      normalization.method,
                      distance.method = "euclidean",
                      KNNGraph_k = 5,
                      data.type = "pseudobulk",
                      sample = NULL,
                      pca = FALSE,
                      ndim = 30,
                      min.samples = 5,
                      min.cells = 10,
                      verbose = TRUE){
   
   mat <- get.matrix(mat,
                     data.type = data.type,
                     ident = ident,
                     sample = sample,
                     min.samples = min.samples,
                     min.cells = min.cells)
   
   # normalized subetted matrix
   norm.mat <- Normalize_data(red.mat@matrix,
                              method = normalization.method)
   nident <- red.mat@ident
   if(pca){
      if(verbose){message("Computing PCA space\n")}
      pr <- custom_prcomp(norm.mat, pca.dims)
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
   
   # set KNNGraph_k
   KNNGraph_k <- min(KNNGraph_k, max(table(nident)))
   # produce KNN
   knn <- compute_KNN(dist,
                      KNNGraph_k = KNNGraph_k)
   
   snn <- compute_snn_graph(dist = dist,
                            KNNGraph_k = KNNGraph_k,
                            knn = knn)
   
   # Get unique cell types
   cell_types <- unique(nident)
   
   # Initialize a list to store proportions for each cell
   cell_proportions <- vector("list", length(nident))
   
   prop.ident.table <- table(nident) |> prop.table()
   
   for (i in seq_along(nident)) {
      # Get the neighbors of the current cell
      neighbors <- knn[i, ]
      
      # Count occurrences of each cell type among neighbors
      neighbor_types <- nident[neighbors]
      type_counts <- table(factor(neighbor_types, levels = cell_types))
      
      score <- type_counts / KNNGraph_k
      # Convert to proportions
      score <- minmax_norm(score,
                           min_value = prop.ident.table[nident[i]],
                           max_value = 1,
                           inverse = F)
      cell_proportions[[i]] <- score
   }
   
   # Convert list to data frame for easier manipulation
   proportion_df <- do.call(rbind, cell_proportions)
   rownames(proportion_df) <- NULL
   colnames(proportion_df) <- cell_types
   
   # Compute the mean proportion of each cell type within each identity group
   aggregate_scores <- aggregate(proportion_df, by = list(nident), FUN = mean)
   colnames(aggregate_scores)[1] <- "celltype"  # Rename grouping column
   
   return(aggregate_scores)
}


dx.singleR.score <- function(pred1, true1,
                             pred2, true2){
   
   pred1 <- singleR.score.tidy(pred1, true1, filter = F)
   pred2 <- singleR.score.tidy(pred2, true2, filter = F)
   pred <- rbind(pred1, pred2)
   
   return(pred)
}

dx.singleR.match <- function(pred1, pred2,
                             true1, true2){
   
   pred1 <- singleR.match.tidy(pred1, true1)
   pred2 <- singleR.match.tidy(pred2, true2)
   
   pred <- dplyr::inner_join(pred1, pred2, by = "true") |>
      dplyr::group_by(true) |>
      dplyr::mutate(score = score.x * score.y) |>
      dplyr::rename("celltype" = "true")
   
   return(pred)
}


dx.bestHit <- function(mat,
                       ident,
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
      
      # original annotations to compare with predictions
      true1 <- sce.list[[a]]$ident
      true2 <- sce.list[[b]]$ident
      
      results <- BiocParallel::bplapply(method,
                                        BPPARAM = param,
                                        function(meth){
                                           
                                           r <- switch(meth,
                                                       "Mutual.Score" = dx.singleR.score(pred1 = pred1,
                                                                                         pred2 = pred2,
                                                                                         true1 = true1,
                                                                                         true2 = true2),
                                                       "Mutual.Match" = dx.singleR.match(pred1 = pred1,
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
      if(m == "Mutual.Score"){
         res <- res |>
            dplyr::group_by(true, celltype) |>
            dplyr::summarize(score = mean(score, na.rm = T))
      } else if (m == "Mutual.Match"){
         res <- res |>
            dplyr::group_by(celltype) |>
            dplyr::summarize(score = mean(score, na.rm = T))
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
      
      return(res)
   }
   )
   names(resli) <- method
   
   return(resli)
}
