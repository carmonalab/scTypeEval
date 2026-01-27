set_parallel_params <- function(ncores = NULL,
                                bparam = NULL,
                                progressbar = TRUE)
{
   if (is.null(ncores)) {
      ncores <- 1
   }
   
   if (ncores > parallel::detectCores()) {
      ncores <- parallel::detectCores()
      message("Using more cores available in this computer, reducing number of cores to ",
              ncores)
   }
   
   # set parallelization parameters
   if (is.null(bparam)) {
      if (ncores > 1) {
         if(.Platform$OS.type == "windows"){
            param <- BiocParallel::SerialParam(progressbar = progressbar)
            warning("Parallel processing not possible for Windows operating system.\n")
         } else {
            param <- BiocParallel::MulticoreParam(workers =  ncores,
                                                  progressbar = progressbar)
         }
      } else {
         param <- BiocParallel::SerialParam(progressbar = progressbar)
      }
   } else {
      param <- bparam
   }
   return(param)
}


split_matrix <- function(mat,
                         ident = NULL,
                         sample = NULL,
                         bparam = BiocParallel::SerialParam()){
   if (length(sample) != ncol(mat)) {
      stop("The length of 'sample' must match the number of columns in 'mat'")
   }
   
   # Get the unique sample groups
   unique_samples <- unique(sample)
   
   # Loop through each unique sample and subset the matrix
   split_list <- BiocParallel::bplapply(unique_samples,
                                        BPPARAM = bparam,
                                        function(s){
                                           new.mat <- mat[, sample == s, drop = FALSE]
                                           new.ident <- ident[sample == s]
                                           
                                           ret <- new("Mat_ident",
                                                      matrix = new.mat,
                                                      group = factor(),
                                                      ident = factor(new.ident),
                                                      sample = s)
                                           
                                           return(ret)
                                        })
   names(split_list) <- as.character(unique_samples)
   
   return(split_list)
}

rm0 <- function(value){
   value[value < 0] <- 0
   return(value)
}


# Function to scale any value to the range [0, 1]
minmax_norm <- function(value, min_value, max_value, inverse = FALSE) {
   norm_value <- (value - min_value) / (max_value - min_value)
   norm_value <- pmax(0, pmin(1, norm_value))  # Clip values to [0,1]
   
   if (inverse) {
      return(1 - norm_value)
   } else {
      return(norm_value)
   }
}




custom_prcomp <- function(norm.mat,
                          ndim = 30, 
                          verbose = TRUE){
   # if ncol or nrow is below given ndim use this number
   ndim <- min(c((dim(norm.mat)-1), ndim))
   
   if(verbose){message("   > Returning ", ndim, " dimensions for PCA")}
   
   # compute PCA
   pr <- irlba::prcomp_irlba(as.matrix(Matrix::t(norm.mat)),
                             n = ndim,
                             center = TRUE,
                             scale. = FALSE)
   
   # keep colnames
   rownames(pr$x) <- colnames(norm.mat)
   
   return(pr)
}

# function to compute variance captured in PC components

var_PCA <- function(pca_embeddings){
   # Compute variance per PC
   n_samples <- nrow(pca_embeddings)  # Number of samples
   variance_per_pc <- colSums(pca_embeddings^2) / (n_samples - 1)
   
   # Compute proportion of variance explained
   total_variance <- sum(variance_per_pc)
   variance_explained <- variance_per_pc / total_variance
   
   return(variance_explained)
}



sample_variable_length_combinations <- function(elements,
                                                min_k = 2,
                                                max_k = NULL,
                                                num_samples = 30,
                                                max_iter = 1000,
                                                seed = NULL) {
   # No global seeding in package code per Bioconductor guidelines
   
   if(is.null(max_k)){
      max_k <- length(elements) - 1
   }
   
   unique_combinations <- list()
   count <- 0
   iter <- 0
   
   # A simple way to sample k values (could be weighted if desired)
   possible_k_values <- min_k:max_k
   
   prop0 <- rep(1, length(elements))
   names(prop0) <- elements
   
   while (count < num_samples) {
      # 1. Randomly pick a size 'k' for the combination
      current_k <- sample(possible_k_values, 1)
      
      # Ensure k is not larger than elements available
      if (current_k > length(elements)) next
      
      
      # 2. Sample 'current_k' elements
      sampled_elements <- sample(elements,
                                 current_k,
                                 replace = FALSE)
      
      # 3. Sort to make it a canonical representation of a combination
      sorted_combination <- sort(sampled_elements)
      current_combo_string <- paste(sorted_combination, collapse = "-")
      
      if (!current_combo_string %in% names(unique_combinations)) {
         unique_combinations[[current_combo_string]] <- sorted_combination
         count <- count + 1
      }
      
      iter <- iter + 1
      
      if(iter >= max_iter){
         break
      }
   }
   
   return(unique_combinations)
}


.general_filtering <- function(mat, # Mat_ident object
                               black.list = NULL,
                               gene.list = NULL,
                               verbose = TRUE){
   norm.mat <- mat@matrix
   # remove blacked listed genes
   if(!is.null(black.list) && verbose){message("   Filtering out black listed genes... \n")}
   norm.mat <- norm.mat[!rownames(norm.mat) %in% black.list,]
   
   # keep only gene list features
   if(verbose){message("   Filtering gene list... \n")}
   norm.mat <- norm.mat[rownames(norm.mat) %in% gene.list,]
   
   # remove rows or columns with only 0
   mat@matrix <- norm.mat
   if(verbose){message("   Filtering empty rows and cols... \n")}
   mat <- filter_empty(mat)
   return(mat)
}

# Helper function to compute centroids for each cluster
compute_centroids <- function(norm.mat, ident) {
   centroid.list <- tapply(seq_len(ncol(norm.mat)), ident, function(idx) {
      Matrix::rowMeans(norm.mat[, idx, drop = FALSE])
   })
   # Combine the list of centroids into a matrix
   centroids <- do.call(cbind, centroid.list) 
   return(centroids)
}

