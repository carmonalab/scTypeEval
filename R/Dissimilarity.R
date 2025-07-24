
dissimilarity_methods <-
   c("WasserStein" = "single-cell",
     "Pseudobulk:Euclidean" = "pseudobulk",
     "Pseudobulk:Cosine" = "pseudobulk",
     "Pseudobulk:Pearson" = "pseudobulk",
     "BestHit:Match" = "pseudobulk",
     "BestHit:Score" = "pseudobulk"
   )


no_dr_ds <- c("BestHit:Match",
              "BestHit:Score")

quiet_transport <- function(...) {
   suppressWarnings(
      suppressMessages(
         capture.output(
            result <- transport::transport(...),
            file = NULL
         )
      )
   )
   result
}


optimal_transport <- function(A, B,
                              p = 2,
                              squared = FALSE){
   # Ensure required packages are available
   if (!requireNamespace("transport", quietly = TRUE)) {
      stop("Please install the 'transport' package first.")
   }
   # Check inputs: A and B must be matrices or data frames with same number of columns (features)
   n <- nrow(A)
   m <- nrow(B)
   # Compute pairwise ground cost matrix between all points in A and B
   # We use Euclidean distance (L2 norm), raised to the power p
   # This is common for Wasserstein-p distances, especially p = 1 or 2
   cost_matrix <- compute_fast_euclidean_2obj(A,B)^p
   # Assign uniform weights to both point sets A and B
   # This assumes empirical measures where each point has equal weight
   a_weights <- rep(1 / n, n)
   b_weights <- rep(1 / m, m)
   # Compute the optimal transport plan using the "shortsimplex" method
   # - "shortsimplex" is based on the revised simplex algorithm with enhancements for efficiency
   # - Recommended by the package for general discrete distributions (class pp)
   # - Handles moderately large problems well
   # - Faster than full simplex but more robust than heuristics like "auction"
   # - Does not require specific conditions (e.g., p=2 or grid-based input)
   suppressWarnings(
      {
         suppressMessages(
            {
               plan <- quiet_transport(a_weights,
                                       b_weights,
                                       costm = cost_matrix,
                                       method = "shortsimplex")
            })})
   
   # Compute the total transportation cost using the optimal plan
   # Multiply transported mass by cost for each (from, to) pair
   total_cost <- sum(plan$mass * cost_matrix[cbind(plan$from, plan$to)])
   
   # Return the Wasserstein distance:
   # - If squared = TRUE, or if p != 2, return the raw p-root of the total cost
   # - Otherwise (default case: p = 2), return the square root for true L2 Wasserstein
   if (squared || p != 2) {
      return(total_cost^(1 / p))
   } else {
      return(sqrt(total_cost))
   }
}


compute_wasserstein <- function(mat,
                                group,
                                transpose = TRUE,
                                p = 2,
                                squared = FALSE,
                                bparam = BiocParallel::SerialParam(),
                                verbose = TRUE){
   
   # split matrix per cell type && sample
   if(verbose){message("Splitting matrices... \n")}
   mat_list <- BiocParallel::bplapply(levels(group),
                                      BPPARAM = bparam,
                                      function(s){
                                         k <- group == s
                                         m <- mat[, k, drop = FALSE]
                                         if(transpose){
                                            m <- Matrix::t(m)
                                         }
                                         return(m)
                                      })
   n <- length(mat_list)
   combs <- utils::combn(n, 2, simplify = FALSE)
   
   # Define the function to compute each pairwise distance
   if(verbose){message("Computing pairwise WasserStein distance... \n")}
   dist_list <- BiocParallel::bplapply(combs,
                                       BPPARAM = bparam,
                                       function(pair) {
                                          i <- pair[1]
                                          j <- pair[2]
                                          d <- optimal_transport(mat_list[[i]], mat_list[[j]],
                                                                 p = p,
                                                                 squared = squared)
                                          list(i = i, j = j, d = d)
                                       })
   
   # Initialize distance matrix
   dist_mat <- matrix(0, n, n)
   rownames(dist_mat) <- colnames(dist_mat) <- levels(group)
   
   # Fill in symmetric matrix
   for (res in dist_list) {
      i <- res$i
      j <- res$j
      d <- res$d
      dist_mat[i, j] <- dist_mat[j, i] <- d
   }
   
   return(as.dist(dist_mat))
}
