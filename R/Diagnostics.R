hirearchy.helper <- function(mat,
                             ident,
                             normalization.method,
                             distance.method = "euclidean",
                             hirearchy.method = "ward.D2",
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
   
   norm.mat <- Normalize_data(mat@matrix,
                              method = normalization.method)
   if(!pca){
      mat <- mat@matrix
   } else {
      pr <- custom_prcomp(norm.mat, ndim)
      mat <- NULL
      norm.mat <- t(pr$x)
   }
   
   dist <- get.distance(mat = mat,
                        norm.mat = norm.mat,
                        distance.method = distance.method)
   hclust_result <- stats::hclust(dist,
                                  method = hirearchy.method)
   
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
   
   norm.mat <- Normalize_data(mat@matrix,
                              method = normalization.method)
   ident <- mat@ident
   if(!pca){
      mat <- mat@matrix
   } else {
      pr <- custom_prcomp(norm.mat, ndim)
      mat <- NULL
      norm.mat <- t(pr$x)
   }
   
   dist <- get.distance(mat = mat,
                        norm.mat = norm.mat,
                        distance.method = distance.method)

   # set KNNGraph_k
   KNNGraph_k <- min(KNNGraph_k, max(table(ident)))
   # produce KNN
   knn <- compute_KNN(dist,
                      KNNGraph_k = KNNGraph_k)
   
   snn <- compute_snn_graph(dist = dist,
                            KNNGraph_k = KNNGraph_k,
                            knn = knn)
   
   # Get unique cell types
   cell_types <- unique(ident)
   
   # Initialize a list to store proportions for each cell
   cell_proportions <- vector("list", length(ident))

   prop.ident.table <- table(ident) |> prop.table()
   
   for (i in seq_along(ident)) {
      # Get the neighbors of the current cell
      neighbors <- knn[i, ]
      
      # Count occurrences of each cell type among neighbors
      neighbor_types <- ident[neighbors]
      type_counts <- table(factor(neighbor_types, levels = cell_types))
      
      score <- type_counts / KNNGraph_k
      # Convert to proportions
      score <- minmax_norm(score,
                           min_value = prop.ident.table[ident[i]],
                           max_value = 1,
                           inverse = F)
      cell_proportions[[i]] <- score
   }
   
   # Convert list to data frame for easier manipulation
   proportion_df <- do.call(rbind, cell_proportions)
   rownames(proportion_df) <- NULL
   colnames(proportion_df) <- cell_types
   
   # Compute the mean proportion of each cell type within each identity group
   aggregate_scores <- aggregate(proportion_df, by = list(ident), FUN = mean)
   colnames(aggregate_scores)[1] <- "celltype"  # Rename grouping column
   
   return(aggregate_scores)
}
