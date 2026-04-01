
distances <- c("euclidean", "cosine", "pearson")

compute_fast_euclidean <- function(norm_mat) {
   # Compute the dot product (squared sums)
   diag_norm <- Matrix::rowSums(norm_mat^2)  # Diagonal terms (squared magnitudes)
   # Compute the squared distance matrix
   dist_squared <- outer(diag_norm, diag_norm, "+") - 2 * Matrix::tcrossprod(norm_mat)
   # Fix numerical precision issues: Negative values close to 0 are set to 0
   dist_squared[dist_squared < 0] <- 0
   # Return the distance matrix
   d <- as.dist(sqrt(dist_squared))
   return(d)
}

compute_fast_euclidean_2obj <- function(a, b) {
   # Squared L2 norms of rows
   a_sq <- Matrix::rowSums(a^2)
   b_sq <- Matrix::rowSums(b^2)
   # Squared Euclidean distances
   dist_sq <- outer(a_sq, b_sq, "+") - 2 * Matrix::tcrossprod(a, b)
   # Fix numerical issues
   dist_sq[dist_sq < 0] <- 0
   # Return distances
   sqrt(dist_sq)
}


compute_cosine <- function(norm_mat) {
   dot_product <- Matrix::tcrossprod(norm_mat)
   magnitude <- sqrt(Matrix::rowSums(norm_mat^2))
   cosine_similarity <- dot_product / outer(magnitude, magnitude)
   as.dist(1 - cosine_similarity)
}

compute_pearson <- function(norm_mat) {
   # Compute row means
   row_means <- Matrix::rowMeans(norm_mat)
   # Center the matrix (subtract row means)
   norm_mat_centered <- norm_mat - row_means
   norm_mat_centered <- Matrix::drop0(norm_mat_centered)  # Keep sparsity
   # Compute covariance matrix
   cov_mat <- Matrix::tcrossprod(norm_mat_centered) / (ncol(norm_mat) - 1)
   # Compute standard deviations
   std_dev <- sqrt(Matrix::diag(cov_mat))
   # Normalize to Pearson correlation matrix
   corr <- cov_mat / (Matrix::tcrossprod(std_dev) + 1e-9)  # Avoid division by zero
   # Convert to distance matrix
   as.dist(1 - corr)
}


# Main function with dispatcher
get_distance <- function(norm_mat = NULL,
                         transpose = TRUE,
                         distance_method = "euclidean",
                         verbose = TRUE) {
   # Handle transposition
   if (transpose) {
      norm_mat <- Matrix::t(norm_mat)
   }
   
   distance_method <- tolower(distance_method)
   
   if(verbose){message("   Running distance for " , distance_method, "... \n")}
   # Dispatcher to call the appropriate function
   dist <- switch(
      distance_method,
      "euclidean" = compute_fast_euclidean(norm_mat),
      "cosine" = compute_cosine(norm_mat),
      "pearson" = compute_pearson(norm_mat),
      stop(distance_method, " is not a supported distance method.")
   )
   
   return(dist)
}
