# Helper function to compute centroids for each cluster
compute_centroids <- function(norm.mat, ident) {
   centroid.list <- tapply(seq_len(ncol(norm.mat)), ident, function(idx) {
      Matrix::rowMeans(norm.mat[, idx, drop = FALSE])
   })
   # Combine the list of centroids into a matrix
   centroids <- do.call(cbind, centroid.list) 
   return(centroids)
}


# min distance between centroids
min.distance_centroids <- function(centroids){
   # Compute pairwise distance matrix between centroids
   centroid_distances <- as.matrix(dist(t(centroids)))
   # Set the diagonal (distance between centroid and itself) to Inf
   diag(centroid_distances) <- Inf
   # Find the minimum distance between distinct centroids
   min_dist_between_centroids <- min(centroid_distances)^2
   
   return(min_dist_between_centroids)
}

# Silhouette Calculation
compute_silhouette <- function(dist, ident) {
   silh <- cluster::silhouette(as.numeric(factor(ident)), dist)
   silh_df <- data.frame(silh)
   silh_df$celltype <- ident
   
   # Calculate the mean silhouette width per cell type
   mean_per_celltype <- silh_df %>%
      dplyr::group_by(celltype) %>%
      dplyr::summarise(mean_sil_width = mean(sil_width, na.rm = TRUE))
   
   return(mean_per_celltype[["mean_sil_width"]])
}

# Modularity Calculation
compute_modularity <- function(dist, ident, KNNGraph_k) {
   knn <- RANN::nn2(as.matrix(dist), k = KNNGraph_k + 1)$nn.idx
   edges <- cbind(rep(1:nrow(knn), each = KNNGraph_k), as.vector(t(knn[, -1])))
   g <- igraph::graph_from_edgelist(edges, directed = FALSE)
   m <- igraph::modularity(g, membership = as.numeric(factor(ident)))
   return(m)
}

# Ward's Method ARI Calculation
compute_ward <- function(dist, ident, method = "ward.D2") {
   hclust_result <- stats::hclust(dist, method = method)
   clusters <- stats::cutree(hclust_result, k = length(unique(ident)))
   w <- mclust::adjustedRandIndex(clusters, ident)
   return(w)
}

# Helper function to compute inertia (sum of squared distances)
# Inertia Calculation
# the lowest more consistency
compute_inertia <- function(norm.mat,
                            ident,
                            centroids) {
   inertia_per_ident <- numeric(length(unique(ident)))
   names(inertia_per_ident) <- unique(ident)
   
   # Compute inertia for each cluster
   for (cluster in unique(ident)) {
      cluster_cells <- norm.mat[,ident == cluster]  # Subset cells for each cluster
      inertia_per_ident[cluster] <- sum(Matrix::colSums((cluster_cells - centroids[,cluster])^2))
   }
   
   return(inertia_per_ident)
}




# Xie-Beni Index Calculation
#A lower Xie-Beni index indicates better clustering quality, meaning that
# clusters are more distinct and compact.
compute_xie_beni <- function(norm.mat,
                             ident,
                             centroids,
                             inertia) {
   # Initialize a vector to store Xie-Beni scores for each ident
   xie_beni_scores <- numeric(length(unique(ident)))
   names(xie_beni_scores) <- unique(ident)
   
   # Iterate over each ident
   for (cluster in unique(ident)) {
      
      # Compute inertia for the current ident
      inertia.sub <- intertia[cluster]
      
      # Compute minimum distance between the current centroid and all other centroids
      other_centroids <- centroids[, setdiff(unique(ident), cluster), drop = FALSE]
      min_dist <- min(as.matrix(dist(t(other_centroids) - centroids[, cluster])))^2
      
      # Number of cells in the current cluster
      n_cells <- sum(ident == cluster)
      
      # Compute Xie-Beni for this ident
      xie_beni_scores[cluster] <- inertia.sub / (n_cells * min_dist)
   }
   
   return(xie_beni_scores)
}

# S_Dbw Index Calculation
compute_s_dbw <- function(norm.mat,
                          ident,
                          centroids) {
   
   # Dispersion
   dispersion <- mean(sapply(unique(ident), function(cluster) {
      cluster_cells <- norm.mat[, ident == cluster, drop = FALSE]  # Subset cells for each cluster
      mean(sqrt(colSums((cluster_cells - centroids[, cluster])^2)))
   }))
   
   # Separation
   separation <- mean(as.dist(1 - cor(t(centroids), method = "pearson")))
   
   s_dbw <- dispersion + separation
   return(s_dbw)
}

# I Index Calculation
compute_i_index <- function(norm.mat, ident, centroids) {
   
   total_centroid <- rowMeans(norm.mat)
   separation <- sum(colSums((centroids - total_centroid)^2))
   cohesion <- sum(sapply(unique(ident), function(cluster) {
      cluster_cells <- colnames(norm.mat)[ident == cluster]  # Subset cells for each cluster
      sum(colSums((norm.mat[, cluster_cells] - centroids[, cluster])^2))
   }))
   
   i_index <- separation / cohesion
   return(i_index)
}


# Main dispatcher for internal validation metrics
calculate_IntVal_metric <- function(norm.mat, 
                                    metric, 
                                    ident, 
                                    dist = NULL, 
                                    centroids = NULL, 
                                    inertia = NULL, 
                                    KNNGraph_k = NULL, 
                                    hclust.method = "ward.D2") {
   
   # Run the requested metric
   switch(metric,
          "silhouette" = compute_silhouette(dist, ident),
          "modularity" = compute_modularity(dist, ident, KNNGraph_k),
          "ward" = compute_ward(dist, ident, hclust.method),
          "inertia" = inertia, 
          "Xie-Beni" = compute_xie_beni(norm.mat, ident, inertia, centroids),
          "S_Dbw" = compute_s_dbw(norm.mat, ident, centroids),
          "I" = compute_i_index(norm.mat, ident, centroids),
          stop(metric, " is not a supported consistency method.")
   )
}
