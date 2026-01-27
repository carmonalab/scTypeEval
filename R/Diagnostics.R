nn.helper <- function(dist,
                      ident,
                      KNNGraph_k,
                      normalize = FALSE
                      ){
   # set KNNGraph_k
   KNNGraph_k <- min(KNNGraph_k, max(table(ident)))
   # produce KNN
   knn <- compute_KNN(dist,
                      KNNGraph_k = KNNGraph_k)
   
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
      if(normalize){
         score <- minmax_norm(score,
                              min_value = prop.ident.table[ident[i]],
                              max_value = 1,
                              inverse = FALSE)
      }
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

