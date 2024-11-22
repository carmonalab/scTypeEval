# Helper functions for each metric
# ident must be a vector with the cell types, names with the cell or sample id!!

compute_silhouette <- function(dist,
                               ident) {
   silh <- cluster::silhouette(as.numeric(factor(ident)),
                               dist)
   silh_df <- data.frame(silh)
   silh_df$celltype <- ident
   
   # Calculate the mean silhouette width per cell type
   mean_per_celltype <- silh_df |>
      dplyr::group_by(celltype) |>
      dplyr::summarise(mean_sil_width = mean(sil_width, na.rm = TRUE))
   
   return(mean_per_celltype)
}

compute_modularity <- function(dist,
                               ident,
                               KNNGraph_k) {
   knn <- RANN::nn2(as.matrix(dist), k = KNNGraph_k + 1)$nn.idx
   edges <- cbind(rep(1:nrow(knn), each = KNNGraph_k), as.vector(t(knn[, -1])))
   g <- igraph::graph_from_edgelist(edges, directed = FALSE)
   igraph::modularity(g, membership = as.numeric(factor(ident)))
}

compute_ward <- function(dist,
                         ident,
                         method = "ward.D2") {
   hclust_result <- stats::hclust(dist,
                                  method = method)
   clusters <- stats::cutree(hclust_result, k = length(unique(ident)))
   w <- mclust::adjustedRandIndex(clusters, ident)
   
   ward_df <- data.frame(ident = unique(ident),
                         ARI = w)
      
   return(ward_df)
}


compute_inertia <- function(norm.mat,
                            ident) {
   
   # Initialize a named vector to store inertia for each ident
   inertia_per_ident <- numeric(length(unique(ident)))
   
   # Compute centroids and inertia for each cluster (cell type)
   for (cluster in unique(ident)) {
      # Subset the cells (columns) belonging to the current ident
      ids <- names(idents)[idents == cluster]
      cluster_points <- norm.mat[, ids, drop = FALSE]
      # Calculate the centroid for the cluster
      centroid <- rowMeans(cluster_points)
      # Calculate within-cluster inertia (sum of squared distances)
      inertia_per_ident[cluster] <- sum(colSums((cluster_points - centroid)^2))
   }
   
   intertia_df <- data.frame(ident = unique(ident),
                             inertia = unlist(inertia_per_ident))
   
   return(inertia_df)
}




calculate_IntVal_metric <- function(dist,
                                    norm.mat,
                                    metric,
                                    ident,
                                    KNNGraph_k,
                                    hclust.method) {
   switch(metric,
          "silhouette" = compute_silhouette(dist, ident),
          "modularity" = compute_modularity(dist, ident, KNNGraph_k),
          "ward" = compute_ward(dist, ident, hclust.method),
          "inertia" = compute_inertia(norm.mat, ident),
          stop(metric_name, " is not a supported consistency method.")
   )
}
