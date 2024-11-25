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

# Modularity Calculation (modified)
# What it measures:
# KNN consistency evaluates the local clustering quality by checking,
# for each cell, how many of its nearest neighbors belong to the same cluster (ident).

compute_modularity <- function(dist,
                               ident,
                               KNNGraph_k = 5) {
   # Create the KNN matrix
   knn <- RANN::nn2(as.matrix(dist), k = KNNGraph_k + 1)$nn.idx
   knn <- knn[, -1]  # Remove self-neighbor
   
   # Initialize a vector to store consistency scores for each cell
   modularities <- numeric(length(ident))
   
   for (i in seq_along(ident)) {
      # Count how many neighbors share the same ident as the current cell
      neighbors <- knn[i, ]
      modularities[i] <- sum(ident[neighbors] == ident[i])
   }
   
   # Normalize by KNNGraph_k to get a proportion
   modularities <- modularities / KNNGraph_k
   
   aggregate_scores <- tapply(modularities, ident, mean)
   
   return(aggregate_scores)
}

# Ward's Method ARI Calculation
compute_ward <- function(dist,
                         ident,
                         method = "ward.D2") {
   
   hclust_result <- stats::hclust(dist,
                                  method = method)
   clusters <- stats::cutree(hclust_result, k = length(unique(ident)))
   # Total cells
   total_cells <- length(ident)
   
   # Proportion of each ident globally
   global_proportion <- table(ident) / total_cells
   
   # Unique cell types
   idents <- unique(ident)
   ward_scores <- numeric(length(idents))
   names(ward_scores) <- idents
   
   for (cluster in idents) {
      # Subset cells belonging to the current ident
      subset_cells <- ident == cluster
      subset_clusters <- clusters[subset_cells]
      
      # Count occurrences of each cluster within this ident
      cluster_counts <- table(subset_clusters)
      
      # Observed proportion in the dominant cluster
      ward_scores[cluster] <- max(cluster_counts) / sum(cluster_counts)
   }
   
   return(ward_scores)
}

# Helper function to compute inertia (sum of squared distances)
# Inertia Calculation
# the lowest more consistency
compute_inertia <- function(norm.mat,
                            ident,
                            centroids = NULL) {
   
   # compute centroids if not provided
   if(is.null(centroids)){
      centroids <- compute_centroids(norm.mat = norm.mat,
                                     ident = ident)
   }
   
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
# it is computed dividing the inertia by the minimum distance to the closer cluster
compute_xie_beni <- function(norm.mat,
                             ident,
                             centroids = NULL,
                             inertia = NULL) {
   
   # compute inertia if not provided
   if(is.null(inertia)){
      intertia <- compute_inertia(norm.mat = norm.mat,
                                  ident = ident,
                                  centroids = centroids)
   }
   
   # Initialize a vector to store Xie-Beni scores for each ident
   xie_beni_scores <- numeric(length(unique(ident)))
   names(xie_beni_scores) <- unique(ident)
   
   # Iterate over each ident
   for (cluster in unique(ident)) {
      # Obtain inertia for the current ident
      inertia.sub <- intertia[cluster]
      
      # Compute minimum distance between the current centroid and all other centroids
      # using Euclidean by default
      other_centroids <- centroids[, setdiff(unique(ident), cluster), drop = FALSE]
      min_dist <- min(dist(t(other_centroids) - centroids[, cluster], method = "euclidean"))^2
      
      # Number of cells in the current cluster
      n_cells <- sum(ident == cluster)
      
      # Compute Xie-Beni for this ident
      xie_beni_scores[cluster] <- inertia.sub / (n_cells * min_dist)
   }
   
   return(xie_beni_scores)
}

# S_Dbw Index Calculation
# 1- Dispersion: Compute the mean Euclidean distance between each cell and the cluster
#  centroid (measuring how "spread out" the cluster is).
# 2 - Separation: Computes the average pairwise dissimilarity between cluster
#  centroids using the Pearson correlation coefficient
# Again the lower the better
compute_s_dbw <- function(norm.mat,
                          ident,
                          centroids = NULL) {
   # compute centroids if not provided
   if(is.null(centroids)){
      centroids <- compute_centroids(norm.mat = norm.mat,
                                     ident = ident)
   }
   
   # Initialize a vector to store S_Dbw values per cluster
   s_dbw_per_cluster <- numeric(length(unique(ident)))
   names(s_dbw_per_cluster) <- unique(ident)
   
   for (cluster in unique(ident)) {
      # Dispersion for the current cluster
      cluster_cells <- norm.mat[, ident == cluster, drop = FALSE]
      dispersion <- mean(sqrt(colSums((cluster_cells - centroids[, cluster])^2)))
      
      # Separation relative to each other clusters
      other_centroids <- centroids[, setdiff(colnames(centroids), cluster), drop = FALSE]
      separation <- mean(1 - apply(other_centroids, 2, 
                                   function(other_centroid){ 
                                      cor(centroids[, cluster], 
                                          other_centroid, 
                                          method = "pearson")
                                   })
      )
      
      
      # Compute S_Dbw for this cluster
      s_dbw_per_cluster[cluster] <- dispersion + separation
   }
   
   return(s_dbw_per_cluster)
}

# I Index Calculation
# The I index is another clustering validation metric designed to evaluate the
# trade-off between separation and cohesion in clusters. It is based on the ratio
# of inter-cluster separation to intra-cluster cohesion:
#    
# 1- Separation measures the dispersion of cluster centroids from the overall dataset centroid.
# 2- Cohesion quantifies the tightness or compactness of each cluster by summing
#  the squared distances between each point and its cluster centroid.

# A higher I index indicates better clustering, as it suggests that clusters are well-separated

compute_i_index <- function(norm.mat,
                            ident,
                            centroids = NULL) {
   # Compute centroids if not provided
   if (is.null(centroids)) {
      centroids <- compute_centroids(norm.mat = norm.mat,
                                     ident = ident)
   }
   
   # Calculate the total dataset centroid
   total_centroid <- rowMeans(norm.mat)
   
   # Initialize a vector to store the I index for each ident
   i_index_scores <- numeric(length(unique(ident)))
   names(i_index_scores) <- unique(ident)
   
   # Iterate over each ident
   for (cluster in unique(ident)) {
      # Compute separation for the current cluster
      separation <- sum((centroids[, cluster] - total_centroid)^2)
      
      # Compute cohesion for the current cluster
      cohesion <- sum(colSums((norm.mat[, ident == cluster] - centroids[, cluster])^2))
      
      # Compute the I index for the current cluster
      i_index_scores[cluster] <- separation / cohesion
   }
   
   return(i_index_scores)
}



# Main dispatcher for internal validation metrics
# Return a vector as long as the number of ident with the consistency value
calculate_IntVal_metric <- function(norm.mat, 
                                    metric, 
                                    ident, 
                                    dist = NULL, 
                                    centroids = NULL, 
                                    inertia = NULL, 
                                    KNNGraph_k = 5, 
                                    hclust.method = "ward.D2") {
   
   # Run the requested metric
   switch(metric,
          "silhouette" = compute_silhouette(dist, ident),
          "modularity" = compute_modularity(dist, ident, KNNGraph_k),
          "ward" = compute_ward(dist, ident, hclust.method),
          "inertia" = if(is.null(inertia)){compute_inertia(norm.mat, ident, centroids)} else {inertia}, 
          "Xie-Beni" = compute_xie_beni(norm.mat, ident, inertia, centroids),
          "S_Dbw" = compute_s_dbw(norm.mat, ident, centroids),
          "I" = compute_i_index(norm.mat, ident, centroids),
          stop(metric, " is not a supported consistency method.")
   )
}
