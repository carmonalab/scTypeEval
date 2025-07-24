
distances <- c("euclidean", "cosine", "pearson")

compute_fast_euclidean <- function(norm.mat) {
   # Compute the dot product (squared sums)
   diag_norm <- Matrix::rowSums(norm.mat^2)  # Diagonal terms (squared magnitudes)
   # Compute the squared distance matrix
   dist_squared <- outer(diag_norm, diag_norm, "+") - 2 * Matrix::tcrossprod(norm.mat)
   # Fix numerical precision issues: Negative values close to 0 are set to 0
   dist_squared[dist_squared < 0] <- 0
   # Return the distance matrix
   d <- as.dist(sqrt(dist_squared))
   return(d)
}

compute_fast_euclidean_2obj <- function(A, B) {
   # Squared L2 norms of rows
   A_sq <- Matrix::rowSums(A^2)
   B_sq <- Matrix::rowSums(B^2)
   # Squared Euclidean distances
   dist_sq <- outer(A_sq, B_sq, "+") - 2 * Matrix::tcrossprod(A, B)
   # Fix numerical issues
   dist_sq[dist_sq < 0] <- 0
   # Return distances
   sqrt(dist_sq)
}


compute_cosine <- function(norm.mat) {
   dot_product <- Matrix::tcrossprod(norm.mat)
   magnitude <- sqrt(Matrix::rowSums(norm.mat^2))
   cosine_similarity <- dot_product / outer(magnitude, magnitude)
   as.dist(1 - cosine_similarity)
}

compute_pearson <- function(norm.mat) {
   # Compute row means
   row_means <- Matrix::rowMeans(norm.mat)
   # Center the matrix (subtract row means)
   norm.mat_centered <- norm.mat - row_means
   norm.mat_centered <- Matrix::drop0(norm.mat_centered)  # Keep sparsity
   # Compute covariance matrix
   cov_mat <- Matrix::tcrossprod(norm.mat_centered) / (ncol(norm.mat) - 1)
   # Compute standard deviations
   std_dev <- sqrt(Matrix::diag(cov_mat))
   # Normalize to Pearson correlation matrix
   corr <- cov_mat / (Matrix::tcrossprod(std_dev) + 1e-9)  # Avoid division by zero
   # Convert to distance matrix
   as.dist(1 - corr)
}


# Main function with dispatcher
get.distance <- function(norm.mat = NULL,
                         transpose = TRUE,
                         distance.method = "euclidean",
                         verbose = TRUE) {
   # Handle transposition
   if (transpose) {
      norm.mat <- Matrix::t(norm.mat)
   }
   
   distance.method <- tolower(distance.method)
   
   if(verbose){message("   Running distance for " , distance.method, "... \n")}
   # Dispatcher to call the appropriate function
   dist <- switch(
      distance.method,
      "euclidean" = compute_fast_euclidean(norm.mat),
      "cosine" = compute_cosine(norm.mat),
      "pearson" = compute_pearson(norm.mat),
      stop(distance.method, " is not a supported distance method.")
   )
   
   return(dist)
}
