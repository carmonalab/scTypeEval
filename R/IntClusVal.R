IntVal_metric <- c("silhouette", "NeighborhoodPurity", "ward.PropMatch",
                  "inertia", "Xie-Beni", "S_Dbw", "I", "Leiden.PropMatch",
                  "modularity")

dist.need <- c("silhouette", "NeighborhoodPurity", "ward.PropMatch",
               "modularity", "Leiden.PropMatch")

knn.need <- c( "NeighborhoodPurity", "modularity", "Leiden.PropMatch")


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
   mean_per_celltype <- silh_df |>
      dplyr::group_by(celltype) |>
      dplyr::summarise(mean_sil_width = mean(sil_width, na.rm = TRUE)) |>
      dplyr::pull(mean_sil_width, name = celltype)
   
   return(mean_per_celltype)
}

compute_KNN <- function(dist, KNNGraph_k){
   # Compute KNN
   knn <- RANN::nn2(as.matrix(dist), k = KNNGraph_k + 1)$nn.idx
   knn <- knn[, -1]  # Remove self-neighbor
   return(knn)
}

# moudularity
compute_snn_graph <- function(dist, KNNGraph_k = 5, knn = NULL) {
   
   if(is.null(knn)){
      knn <- compute_KNN(dist = dist, KNNGraph_k = KNNGraph_k)
   }
   
   # Initialize adjacency matrix
   n <- nrow(as.matrix(dist))
   adj_matrix <- matrix(0, n, n)
   
   # Count shared neighbors
   for (i in seq_len(n)) {
      for (j in knn[i, ]) {
         shared_neighbors <- length(intersect(knn[i, ], knn[j, ]))
         adj_matrix[i, j] <- shared_neighbors
         adj_matrix[j, i] <- shared_neighbors  # Ensure symmetry
      }
   }
   
   # Create graph object
   g <- igraph::graph_from_adjacency_matrix(adj_matrix,
                                            mode = "undirected",
                                            weighted = TRUE,
                                            diag = FALSE)
   return(g)
}


compute_modularity_global <- function(dist,
                                      ident,
                                      KNNGraph_k = 5,
                                      knn = NULL) {
   # Create a graph object
   g <- compute_snn_graph(dist = dist,
                          KNNGraph_k = KNNGraph_k,
                          knn = knn)
   
   # Compute modularity
   modularity_score <- igraph::modularity(g, membership = as.numeric(factor(ident)))
   names(modularity_score) <- "global"
   
   return(modularity_score)
}

# Uses the modularity matrix to calculate contributions for specific cell types directly.

# This code computes the modularity contribution per cell type using a graph-based approach,
# where the input distance matrix is converted into a graph, and community detection is
# approximated based on a predefined identity (ident) assigned to each node.

compute_modularity_pct <- function(dist,
                                   ident,
                                   KNNGraph_k = 5,
                                   knn = NULL) {
   # Create a graph object
   g <- compute_snn_graph(dist = dist,
                          KNNGraph_k = KNNGraph_k,
                          knn = knn)
   
   # Get the modularity matrix
   mod_matrix <- igraph::modularity_matrix(g)
   
   # Compute modularity contributions per cell type
   cell_types <- unique(ident)
   modularity_per_type <- sapply(cell_types, function(ct) {
      type_cells <- which(ident == ct)  # Indices of cells in this type
      num_type_cells <- length(type_cells)  # Number of cells in this type
      
      # Modularity contribution normalized by the number of type cells
      type_modularity <- (sum(mod_matrix[type_cells, type_cells]) / num_type_cells) /
                           (2 * igraph::ecount(g))
      return(type_modularity)
   })
   
   # Return a named vector of modularity contributions
   names(modularity_per_type) <- cell_types
   return(modularity_per_type)
}


# NeighborhoodPurity 
# What it measures:
# KNN consistency evaluates the local clustering quality by checking,
# for each cell, how many of its nearest neighbors belong to the same cluster (ident).

compute_NeighborhoodPurity <- function(dist,
                               ident,
                               KNNGraph_k = 5,
                               knn = NULL) {
   if(is.null(knn)){
      knn <- compute_KNN(dist = dist, KNNGraph_k = KNNGraph_k)
   }
   
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
   aggregate_scores <- aggregate_scores[!is.na(aggregate_scores)]
   
   # convert to numeric vector
   an <- names(aggregate_scores)
   aggregate_scores <- as.vector(aggregate_scores)
   names(aggregate_scores) <- an
   
   return(aggregate_scores)
}


custom_PropMatch <- function(ident, clusters){
   # Unique cell types
   idents <- unique(ident)
   scores <- numeric(length(idents))
   names(scores) <- idents
   
   
   prop.ident <- table(ident) 
   prop.ident.table <- prop.ident |> prop.table()
   
   for (ct in idents) {
      # Subset cells belonging to the current ident
      subset_cells <- ident == ct
      subset_clusters <- clusters[subset_cells]
      
      # Count occurrences of each cluster within this ident
      cluster_counts <- table(subset_clusters)
      major_cluster <- names(which.max(cluster_counts))
      
      # Observed proportion in the dominant cluster and normalize for the expected proportion
      score <- (cluster_counts[major_cluster] / prop.ident[ct]) 
      
      score <- minmax_norm(score,
                           min_value = prop.ident.table[ct],
                           max_value = 1,
                           inverse = F)
      scores[ct] <- score
   }
   
   return(scores)
}


# Ward's Method ARI Calculation
compute_ward <- function(dist,
                         ident,
                         method = "ward.D2") {
   
   hclust_result <- stats::hclust(dist,
                                  method = method)
   clusters <- stats::cutree(hclust_result, k = length(unique(ident)))
   scores <- custom_PropMatch(ident, clusters)
   return(scores)
}


compute_leiden <- function(dist,
                           ident,
                           KNNGraph_k = 5,
                           knn = NULL,
                           initial_resolution = 1e-3,
                           max_iter = 500) {
   
   target_clusters <- nlevels(ident)

   g <- compute_snn_graph(dist = dist,
                          KNNGraph_k = KNNGraph_k,
                          knn = knn)

   # Adjust resolution to match the target number of clusters
   resolution <- initial_resolution
   for (i in seq_len(max_iter)) {
      leiden <- igraph::cluster_leiden(g, resolution = resolution)
      clusters <- igraph::membership(leiden)
      
      num_clusters <- length(unique(clusters))
      if (num_clusters == target_clusters) {
         break
      }
      
      # Dynamically adjust resolution (simple heuristic)
      if (num_clusters != target_clusters) {
         resolution <- resolution * (num_clusters/target_clusters)
      } 
   }
   
   if (num_clusters != target_clusters){
      stop("Leiden could not find expected numbrer of clusters, change KNNGraph_k")
   }
   
   # Step 4: Compute probability match
   scores <- custom_PropMatch(ident, clusters)
   return(scores)
}


# Helper function to compute inertia (sum of squared distances)
# Inertia Calculation
# the lowest more consistency
# Global: summing the per-cluster inertia values gives the overall inertia for the dataset.
# The mean of per-cluster inertia values would divide by the number of clusters,
# but the global inertia is the sum across all clusters, not the mean.
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
# Global: summing the per-cluster Xie-Beni values does not yield the correct global Xie-Beni index.
# The denominator of the global Xie-Beni includes the minimum distance between any two clusters
# and the total number of cells. Summing per-cluster Xie-Beni values will distort the result
# because the minimum distance is not shared across clusters
# Averaging per-cluster Xie-Beni values would still not reflect the global denominator 

compute_xie_beni <- function(norm.mat,
                             ident,
                             centroids = NULL,
                             inertia = NULL) {
   
   # compute inertia if not provided
   if(is.null(inertia)){
      inertia <- compute_inertia(norm.mat = norm.mat,
                                  ident = ident,
                                  centroids = centroids)
   }
   
   # Initialize a vector to store Xie-Beni scores for each ident
   xie_beni_scores <- numeric(length(unique(ident)))
   names(xie_beni_scores) <- unique(ident)
   
   # Iterate over each ident
   for (cluster in unique(ident)) {
      # Obtain inertia for the current ident
      inertia.sub <- inertia[cluster]
      
      # Compute minimum distance between the current centroid and all other centroids
      # using Euclidean by default
      # Check if there are more than one cluster to compute minimum distance
      if (length(unique(ident)) > 2) {
         # Compute minimum distance between the current centroid and all other centroids
         # using Euclidean distance by default
         other_centroids <- centroids[, setdiff(unique(ident), cluster), drop = FALSE]
         min_dist <- min(dist(t(other_centroids) - centroids[, cluster], method = "euclidean"))^2
      } else {
         # If there are only two clusters, compute the distance directly between the two centroids
         # Compute distance between the centroids (squared Euclidean distance)
         min_dist <- sum((centroids[, setdiff(unique(ident), cluster)] - centroids[, cluster])^2)
      }
      
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

# Global: summing per-cluster S_Dbw values will not yield the global S_Dbw index.
# The global dispersion and global separation must be averaged across clusters,
# rather than summed. Summing will overestimate the result and misrepresent the
# balance between dispersion and separation.

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
      dispersion <- mean(sqrt(Matrix::colSums((cluster_cells - centroids[, cluster])^2)))
      
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

# Global: summing per-cluster I index values does not yield the global I index.
# The global I index is a ratio of total separation to total cohesion.
# Simply summing the per-cluster ratios will not preserve this relationship.

compute_i_index <- function(norm.mat,
                            ident,
                            centroids = NULL) {
   # Compute centroids if not provided
   if (is.null(centroids)) {
      centroids <- compute_centroids(norm.mat = norm.mat,
                                     ident = ident)
   }
   
   # Calculate the total dataset centroid
   total_centroid <- Matrix::rowMeans(norm.mat)
   
   # Initialize a vector to store the I index for each ident
   i_index_scores <- numeric(length(unique(ident)))
   names(i_index_scores) <- unique(ident)
   
   # Iterate over each ident
   for (cluster in unique(ident)) {
      # Compute separation for the current cluster
      separation <- sum((centroids[, cluster] - total_centroid)^2)
      
      # Compute cohesion for the current cluster
      cohesion <- sum(Matrix::colSums((norm.mat[, ident == cluster] - centroids[, cluster])^2))
      
      # Compute the I index for the current cluster
      i_index_scores[cluster] <- separation / cohesion
   }
   
   return(i_index_scores)
}



# Main dispatcher for internal validation metrics
# Return a vector as long as the number of ident with the consistency value
calculate_IntVal_metric <- function(mat = NULL,
                                    norm.mat = NULL,
                                    metrics,  # Accept a vector of metric names
                                    dist.method = "euclidean",
                                    ident, 
                                    dist = NULL, 
                                    centroids = NULL, 
                                    inertia = NULL, 
                                    KNNGraph_k = 5, 
                                    hclust.method = "ward.D2",
                                    verbose = T) {
   
   # Validate metrics input
   if (!all(metrics %in% IntVal_metric)) {
      stop("One or more requested metrics are not supported.")
   }
   
   # Precompute distance if needed
   
   if (is.null(dist) && any(metrics %in% dist.need)) {
      if(verbose){message("Computing distances...\n")}
      dist <- get.distance(mat = mat,
                           norm.mat = norm.mat,
                           dist.method = dist.method)
   }
   
   if(any(metrics %in% knn.need)){
      # set KNNGraph_k
      knn <- compute_KNN(dist, KNNGraph_k = KNNGraph_k)
   }
   
   if(is.null(centroids) && any(!metrics %in% dist.need)){
      if(verbose){message("Computing centroids...\n")}
      centroids <- compute_centroids(norm.mat, ident)
   }
   
   
   # Precompute inertia if needed
   if ("inertia" %in% metrics && is.null(inertia)) {
      if(verbose){message("Computing inertia...\n")}
      inertia <- compute_inertia(norm.mat = norm.mat,
                                 ident = ident,
                                 centroids = centroids)
   }
   
   # Initialize results list
   results <- list()
   
   # Run requested metrics
   for (metric in metrics) {
      results[[metric]] <- switch(metric,
                                  "silhouette" = compute_silhouette(dist, ident),
                                  "NeighborhoodPurity" = compute_NeighborhoodPurity(dist, ident, KNNGraph_k, knn),
                                  "modularity" = compute_modularity_global(dist, ident, KNNGraph_k, knn),
                                  "ward.PropMatch" = compute_ward(dist, ident, hclust.method),
                                  "Leiden.PropMatch" = compute_leiden(dist, ident, KNNGraph_k, knn),
                                  "inertia" = inertia, 
                                  "Xie-Beni" = compute_xie_beni(norm.mat, ident, centroids, inertia),
                                  "S_Dbw" = compute_s_dbw(norm.mat, ident, centroids),
                                  "I" = compute_i_index(norm.mat, ident, centroids),
                                  stop(metric, " is not a supported consistency method.")
      )
   }
   
   return(results)
}

