distance_methods <- c(
   "euclidean",
   "maximum",
   "manhattan",
   "canberra",
   "binary",
   "minkowski",
   "Jaccard",
   "Weighted_Jaccard",
   "gower",
   "bray-curtis",
   "cosine",
   "pearson"
)

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





# Define sub-functions for specific distance methods
compute_stat_dist <- function(norm.mat, dist.method) {
   stats::dist(norm.mat, method = dist.method)
}

compute_jaccard <- function(mat) {
   mat[mat != 0] <- 1  # Convert to binary
   as.dist(1 - Matrix::tcrossprod(mat) / Matrix::rowSums(mat + t(mat)))  # Jaccard distance calculation without extra dependencies
}


# Function to compute the Weighted Jaccard distance
weighted_jaccard <- function(x, y) {
   sum_pmin <- sum(pmin(x, y))
   sum_pmax <- sum(pmax(x, y))
   return(1 - (sum_pmin / sum_pmax))  # 1 - Jaccard for a distance measure
}

compute_weighted_jaccard <- function(norm.mat) {
   n <- nrow(norm.mat)
   dist_matrix <- matrix(0, n, n)
   for (i in seq_len(n)) {
      for (j in seq_len(n)) {
         dist_matrix[i, j] <- weighted_jaccard(norm.mat[i, ], norm.mat[j, ])
      }
   }
   as.dist(dist_matrix)
}

compute_gower <- function(mat, norm.mat) {
   mat[mat != 0] <- 1  # Convert to binary
   cdata <- cbind(mat, norm.mat)
   cluster::daisy(cdata, metric = "gower")  # `cluster` is already used in Bioconductor
}

compute_bray_curtis <- function(norm.mat) {
   row_sums <- Matrix::rowSums(norm.mat)
   as.dist(1 - Matrix::tcrossprod(norm.mat) / outer(row_sums, row_sums, "+"))
}

compute_cosine <- function(norm.mat) {
   dot_product <- Matrix::tcrossprod(norm.mat)
   magnitude <- sqrt(Matrix::rowSums(norm.mat^2))
   cosine_similarity <- dot_product / outer(magnitude, magnitude)
   as.dist(1 - cosine_similarity)
}

compute_pearson <- function(norm.mat) {
   corr <- stats::cor(t(norm.mat), method = "pearson")
   as.dist(1 - corr)
}


# Main function with dispatcher
get.distance <- function(mat = NULL,
                         norm.mat = NULL,
                         transpose = TRUE,
                         dist.method = "euclidean") {
   # Handle transposition
   if (transpose) {
      if(!is.null(mat)){mat <- Matrix::t(mat)}
      if(!is.null(norm.mat)){norm.mat <- Matrix::t(norm.mat)}
   }
   
   # Dispatcher to call the appropriate function
   dist <- switch(
      dist.method,
      "euclidean" = compute_fast_euclidean(norm.mat),
      "maximum" = compute_stat_dist(norm.mat, dist.method),
      "manhattan" = compute_stat_dist(norm.mat, dist.method),
      "canberra" = compute_stat_dist(norm.mat, dist.method),
      "binary" = compute_stat_dist(norm.mat, dist.method),
      "minkowski" = compute_stat_dist(norm.mat, dist.method),
      "Jaccard" = compute_jaccard(mat),
      "Weighted_Jaccard" = compute_weighted_jaccard(norm.mat),
      "gower" = compute_gower(mat, norm.mat),
      "bray-curtis" = compute_bray_curtis(norm.mat),
      "cosine" = compute_cosine(norm.mat),
      "pearson" = compute_pearson(norm.mat),
      
      stop(dist.method, " is not a supported distance method.")
   )
   
   return(dist)
}
