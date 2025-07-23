IntVal_metric <- c("silhouette", "NeighborhoodPurity",
                   "ward.PropMatch", "Orbital.medoid",
                   "Average.similarity")

dist.need <- "silhouette|NeighborhoodPurity|ward|Orbital.medoid|Average.similarity"

knn.need <- "NeighborhoodPurity"

# compute medoid based on distances
# expecting symmetric distances
compute_medoid_indices <- function(ident,
                                   dist) {
   dist_mat <- as.matrix(dist)  # full symmetric matrix
   medoid_indices <- tapply(seq_along(ident), ident, function(idx) {
      sub_dists <- dist_mat[idx, idx, drop = FALSE]
      total_dist <- rowSums(sub_dists)
      idx[which.min(total_dist)]
   })
   return(unlist(medoid_indices[!is.na(medoid_indices)]))
}

compute_orbital_proportion <- function(norm.mat = NULL,
                                       ident,
                                       dist,
                                       centroids = NULL,
                                       use = "medoid") {
   use <- use[1]
   
   medoid_idx <- NULL
   d_to_reps <- NULL
   
   if (use == "centroid") {
      # Get or compute centroids
      reps <- if (is.null(centroids)) {
         compute_centroids(norm.mat, ident)
      } else {
         centroids
      }
      
      # Compute pairwise Euclidean distances from all cells to centroids
      x <- as.matrix(Matrix::t(norm.mat))
      y <- as.matrix(reps)
      d_to_reps <- sqrt(
         outer(rowSums(x^2), colSums(y^2), "+") - 2 * x %*% y
      )
      colnames(d_to_reps) <- colnames(y)
      
   } else if (use == "medoid") {
      dist_mat <- as.matrix(dist)
      medoid_idx <- compute_medoid_indices(ident, dist)
      
      # Subset distance matrix: rows = all cells, columns = medoids
      d_to_reps <- dist_mat[, medoid_idx, drop = FALSE]
      colnames(d_to_reps) <- names(medoid_idx)
   }
   
   # Rows: samples, Columns: cluster representatives
   own_clusters <- as.character(ident)
   own_idx <- match(own_clusters, colnames(d_to_reps))
   
   own_dist <- d_to_reps[cbind(seq_along(own_clusters), own_idx)]
   
   min_other_dist <- vapply(seq_len(nrow(d_to_reps)),
                            function(i) {
                               row <- d_to_reps[i, ]
                               own <- own_idx[i]
                               min(row[-own], na.rm = TRUE)
                            },
                            numeric(1))
   
   closer_to_own <- own_dist < min_other_dist
   
   # If using medoids, exclude the medoid samples themselves
   if (use == "medoid") {
      medoid_cells <- as.vector(medoid_idx)
      ident <- ident[-medoid_cells]
      closer_to_own <- closer_to_own[-medoid_cells]
   }
   
   # Compute and return the proportion vector
   prop_vector <- tapply(closer_to_own, ident, mean)
   prop_vector <- setNames(as.numeric(prop_vector), names(prop_vector))
   prop_vector <- prop_vector[!is.na(prop_vector)]
   return(prop_vector)
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
   # Convert distance object to a matrix (if necessary)
   dist_matrix <- as.matrix(dist)
   
   # For each row, find the indices of the k nearest neighbors
   knn <- apply(dist_matrix, 1, function(row) {
      order(row)[1:(KNNGraph_k + 1)]  # Include self-neighbor
   })
   # Transpose the result and remove the first column (self-neighbor)
   knn <- t(knn)[, -1]
   return(knn)
}

# moudularity
compute_snn_graph <- function(dist, KNNGraph_k = 5, knn = NULL) {
   
   if(is.null(knn)){
      knn <- compute_KNN(dist = dist, KNNGraph_k = KNNGraph_k)
   }
   
   # Initialize adjacency matrix
   n <- attr(dist, "Size")
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
                                      knn = NULL,
                                      snn = NULL) {
   # Create a graph object
   if(is.null(snn)){
      snn <- compute_snn_graph(dist = dist,
                               KNNGraph_k = KNNGraph_k,
                               knn = knn)
   }
   
   # Compute modularity
   modularity_score <- igraph::modularity(snn, membership = as.numeric(factor(ident)))
   names(modularity_score) <- "global"
   
   return(modularity_score)
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
      knn <- compute_KNN(dist = dist,
                         KNNGraph_k = KNNGraph_k)
   }
   
   # Normalize by expected proportion
   prop.ident <- table(ident) 
   prop.ident.table <- prop.ident |> prop.table()
   # Initialize a vector to store consistency scores for each cell
   modularities <- numeric(length(ident))
   
   for (i in seq_along(ident)) {
      # Count how many neighbors share the same ident as the current cell
      neighbors <- knn[i, ]
      score <- sum(ident[neighbors] == ident[i]) / KNNGraph_k
      # normalize by random chance
      score <- minmax_norm(score,
                           min_value = prop.ident.table[ident[i]],
                           max_value = 1,
                           inverse = F)
      
      modularities[i] <- score
   }
   
   # compute mean per cell type
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
      # normalize by random chance
      score <- minmax_norm(score,
                           min_value = prop.ident.table[ct],
                           max_value = 1,
                           inverse = F)
      scores[ct] <- score
   }
   
   return(scores)
}

compute_PropMatch <- function(ident,
                              clusters,
                              label.conservation = c("PropMatch", "NMI", "ARI")){
   
   label.conservation <- label.conservation[1]
   scores <- switch(label.conservation,
                    "PropMatch" = custom_PropMatch(ident,
                                                   clusters),
                    "NMI" = aricode::NMI(factor(ident),
                                         factor(clusters)),
                    "ARI" = aricode::ARI(factor(ident),
                                         factor(clusters)),
                    stop(label.conservation, " is not a supported label conservation metrics.")
   )
   # change name if global score
   if(label.conservation %in% c("NMI", "ARI")){
      names(scores) <- "global"
   }
   return(scores)
}

# proportion matching
compute_ward <- function(dist,
                         ident,
                         method = "ward.D2",
                         label.conservation = "PropMatch") {
   
   hclust_result <- stats::hclust(dist,
                                  method = method)
   clusters <- stats::cutree(hclust_result,
                             k = length(unique(ident)))
   scores <- compute_PropMatch(ident,
                               clusters,
                               label.conservation)
   return(scores)
}




# Main dispatcher for internal validation metrics
# Return a vector as long as the number of ident with the consistency value
calculate_IntVal_metric <- function(mat = NULL,
                                    norm.mat = NULL,
                                    metrics,  # Accept a vector of metric names
                                    distance.method = "euclidean",
                                    ident, 
                                    dist = NULL, 
                                    centroids = NULL, 
                                    inertia = NULL, 
                                    knn = NULL,
                                    snn = NULL,
                                    KNNGraph_k = 5, 
                                    hclust.method = "ward.D2",
                                    verbose = T) {
   
   # Validate metrics input
   if (!all(metrics %in% IntVal_metric)) {
      stop("One or more requested metrics are not supported.")
   }
   
   
   # Precompute distance if needed
   
   if (is.null(dist) && any(grepl(dist.need, metrics))) {
      if(verbose){message("Computing distances...\n")}
      dist <- get.distance(mat = mat,
                           norm.mat = norm.mat,
                           distance.method = distance.method)
   }
   
   if(any(grepl(knn.need, metrics))){
      # set KNNGraph_k
      KNNGraph_k <- min(KNNGraph_k, max(table(ident)))
      # produce KNN
      knn <- compute_KNN(dist,
                         KNNGraph_k = KNNGraph_k)
   }
   
   # Initialize results list
   results <- list()
   
   # Run requested metrics
   for (metric in metrics) {
      results[[metric]] <- 
         switch(metric,
                "silhouette" = compute_silhouette(dist,
                                                  ident),
                "NeighborhoodPurity" = compute_NeighborhoodPurity(dist,
                                                                  ident,
                                                                  KNNGraph_k,
                                                                  knn),
                "ward.PropMatch" = compute_ward(dist,
                                                ident,
                                                hclust.method,
                                                label.conservation = "PropMatch"),
                "Orbital.medoid" = compute_orbital_proportion(ident = ident,
                                                              dist = dist,
                                                              use = "medoid"),
                "Average.similarity" = compute_averageSimilarity(),
                stop(metric, " is not a supported consistency method.")
         )
   }
   
   return(results)
}

