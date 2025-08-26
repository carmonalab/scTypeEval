IntVal_metric <- c("silhouette", "NeighborhoodPurity",
                   "ward.PropMatch", "Orbital.medoid",
                   "Average.similarity")
knn.need <- "NeighborhoodPurity"

# compute medoid based on distances
# expecting symmetric distances
compute_medoid_indices <- function(dist,
                                   ident) {
   dist_mat <- as.matrix(dist)  # full symmetric matrix
   medoid_indices <- tapply(seq_along(ident), ident, function(idx) {
      sub_dists <- dist_mat[idx, idx, drop = FALSE]
      total_dist <- rowSums(sub_dists)
      idx[which.min(total_dist)]
   })
   return(unlist(medoid_indices[!is.na(medoid_indices)]))
}

compute_orbital_proportion <- function(dist,
                                       ident,
                                       normalize = FALSE) {
   
   dist_mat <- as.matrix(dist)

   medoid_idx <- compute_medoid_indices(dist, ident)
   
   # Subset distance matrix: rows = all samples, columns = medoids
   d_to_reps <- dist_mat[, medoid_idx, drop = FALSE]
   colnames(d_to_reps) <- names(medoid_idx)
   
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
   
   # Exclude the medoid cells themselves
   keep_idx <- setdiff(seq_along(closer_to_own), medoid_idx)
   ident <- ident[keep_idx]
   closer_to_own <- closer_to_own[keep_idx]
   
   # Compute raw proportion per cluster
   prop_vector <- tapply(closer_to_own, ident, mean)
   prop_vector <- setNames(as.numeric(prop_vector), names(prop_vector))
   prop_vector <- prop_vector[!is.na(prop_vector)]
   
   # Normalization by expected proportion (like compute_NeighborhoodPurity)
   if(normalize) {
      prop.ident <- table(ident)
      prop.ident.table <- prop.ident |> prop.table()
      
      for (cl in names(prop_vector)) {
         expected <- prop.ident.table[cl]
         prop_vector[cl] <- minmax_norm(prop_vector[cl],
                                        min_value = expected,
                                        max_value = 1,
                                        inverse = FALSE)
      }
   }
   
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


# NeighborhoodPurity 
# What it measures:
# KNN consistency evaluates the local clustering quality by checking,
# for each cell, how many of its nearest neighbors belong to the same cluster (ident).

compute_NeighborhoodPurity <- function(dist,
                                       ident,
                                       KNNGraph_k = 5,
                                       knn = NULL,
                                       normalize = FALSE) {
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
      if(normalize){
         score <- minmax_norm(score,
                              min_value = prop.ident.table[ident[i]],
                              max_value = 1,
                              inverse = F)
      }
      
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


custom_PropMatch <- function(ident, clusters, normalize = FALSE){
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
      if(normalize){
         score <- minmax_norm(score,
                              min_value = prop.ident.table[ct],
                              max_value = 1,
                              inverse = F)
      }
      scores[ct] <- score
   }
   
   return(scores)
}

compute_PropMatch <- function(ident,
                              clusters,
                              label.conservation = c("PropMatch", "NMI", "ARI"),
                              normalize = FALSE){
   
   label.conservation <- label.conservation[1]
   scores <- switch(label.conservation,
                    "PropMatch" = custom_PropMatch(ident,
                                                   clusters,
                                                   normalize),
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
                         label.conservation = "PropMatch",
                         normalize = FALSE) {
   
   hclust_result <- stats::hclust(dist,
                                  method = method)
   clusters <- stats::cutree(hclust_result,
                             k = length(unique(ident)))
   scores <- compute_PropMatch(ident,
                               clusters,
                               label.conservation,
                               normalize = normalize)
   return(scores)
}

# similar to one-label silhouette
compute_averageSimilarity <- function(dist, ident) {
   dist_mat <- as.matrix(dist)
   ident <- as.character(ident)  # ensure consistent indexing
   
   unique_idents <- unique(ident)
   ratio_vector <- numeric(length(ident))
   
   for (i in seq_along(ident)) {
      this_ident <- ident[i]
      
      same_group_idx <- which(ident == this_ident & seq_along(ident) != i)
      other_group_idx <- which(ident != this_ident)
      
      # Handle edge cases: if only one cell in group, skip
      if (length(same_group_idx) == 0 || length(other_group_idx) == 0) {
         ratio_vector[i] <- NA
         next
      }
      
      avg_own_dist <- mean(dist_mat[i, same_group_idx])
      avg_other_dist <- mean(dist_mat[i, other_group_idx])
      
      # Compute normalized score: higher is better (self-similar and dissimilar to others)
      ratio_vector[i] <- 1 - (avg_own_dist / (avg_own_dist + avg_other_dist))
   }
   
   # Aggregate: mean per cluster
   similarity_score <- tapply(ratio_vector, ident, mean, na.rm = TRUE)
   
   similarity_score <- similarity_score[!is.na(similarity_score)]
   similarity_score <- setNames(as.numeric(similarity_score), names(similarity_score))
   
   return(similarity_score)
}



# Main dispatcher for internal validation metrics
# Return a vector as long as the number of ident with the consistency value
calculate_IntVal_metric <- function(dist,
                                    ident,
                                    metrics, 
                                    KNNGraph_k = 5, 
                                    hclust.method = "ward.D2",
                                    normalize = FALSE) {
   
   # Validate metrics input
   if (!all(metrics %in% IntVal_metric)) {
      stop("One or more requested metrics are not supported.")
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
                "NeighborhoodPurity" = compute_NeighborhoodPurity(dist = dist,
                                                                  ident = ident,
                                                                  KNNGraph_k = KNNGraph_k,
                                                                  knn = knn,
                                                                  normalize = normalize),
                "ward.PropMatch" = compute_ward(dist,
                                                ident,
                                                hclust.method,
                                                label.conservation = "PropMatch",
                                                normalize = normalize),
                "Orbital.medoid" = compute_orbital_proportion(ident = ident,
                                                              dist = dist,
                                                              normalize = normalize),
                "Average.similarity" = compute_averageSimilarity(dist, ident),
                stop(metric, " is not a supported consistency method.")
         )
   }
   
   return(results)
}

