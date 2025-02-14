distance_methods <- c(
   "euclidean",
   "EMD",
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

#earth mover distance (emd)
compute_emd <- function(norm.mat,
                        dist = "euclidean",
                        bw = NULL,  # Set NULL to auto-calculate per pseudobulk
                        grid_size = 100,
                        max.iter = 1000) {
   n <- nrow(norm.mat)  # Each row is a pseudobulk sample
   
   # Convert the matrix to dense if it's sparse
   if (inherits(norm.mat, "dgCMatrix")) {
      norm.mat <- as.matrix(norm.mat)
   }
   
   # Define a grid for KDE estimation
   min_val <- min(norm.mat, na.rm = TRUE)
   max_val <- max(norm.mat, na.rm = TRUE)
   grid <- seq(min_val, max_val, length.out = grid_size)
   
   # Preallocate KDE matrix
   kde_matrix <- matrix(0, nrow = grid_size, ncol = n)  
   rownames(kde_matrix) <- paste0("KDE_", seq_len(grid_size))
   colnames(kde_matrix) <- rownames(norm.mat)  # Retain original sample names
   
   # Compute KDE for each pseudobulk (each row)
   for (i in 1:n) {
      data_i <- as.numeric(norm.mat[i, ])  # Ensure it's numeric
      h_i <- if (is.null(bw)) stats::bw.nrd0(data_i) else bw  # Compute bandwidth if needed
      kde <- ks::kde(data_i, h = h_i, eval.points = grid)  
      kde_matrix[, i] <- kde$estimate / sum(kde$estimate)  # Normalize
   }
   
   # Preallocate distance vector
   emd_distances <- numeric(n * (n - 1) / 2)
   idx <- 1
   
   # Compute pairwise EMD
   for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
         emd_distances[idx] <- emdist::emd2d(kde_matrix[, i, drop = F],
                                           kde_matrix[, j, drop = F],
                                           dist = dist,
                                           max.iter = max.iter)
         idx <- idx + 1
      }
   }
   
   # Convert to `dist` object with rownames
   dist_object <- structure(
      emd_distances,
      Size = n,
      Labels = rownames(norm.mat),  # Keep original rownames
      class = "dist",
      Diag = FALSE,
      Upper = FALSE
   )
   return(dist_object)
}




# Define sub-functions for specific distance methods
compute_stat_dist <- function(norm.mat, distance.method) {
   stats::dist(norm.mat, method = distance.method)
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
get.distance <- function(mat = NULL,
                         norm.mat = NULL,
                         transpose = TRUE,
                         distance.method = "euclidean") {
   # Handle transposition
   if (transpose) {
      if(!is.null(mat)){mat <- Matrix::t(mat)}
      if(!is.null(norm.mat)){norm.mat <- Matrix::t(norm.mat)}
   }
   
   # Dispatcher to call the appropriate function
   dist <- switch(
      distance.method,
      "euclidean" = compute_fast_euclidean(norm.mat),
      "EMD" = compute_emd(norm.mat),
      "maximum" = compute_stat_dist(norm.mat, distance.method),
      "manhattan" = compute_stat_dist(norm.mat, distance.method),
      "canberra" = compute_stat_dist(norm.mat, distance.method),
      "binary" = compute_stat_dist(norm.mat, distance.method),
      "minkowski" = compute_stat_dist(norm.mat, distance.method),
      "Jaccard" = compute_jaccard(mat),
      "Weighted_Jaccard" = compute_weighted_jaccard(norm.mat),
      "gower" = compute_gower(mat, norm.mat),
      "bray-curtis" = compute_bray_curtis(norm.mat),
      "cosine" = compute_cosine(norm.mat),
      "pearson" = compute_pearson(norm.mat),
      
      stop(distance.method, " is not a supported distance method.")
   )
   
   return(dist)
}
