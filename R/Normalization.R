# we will offer 3 types of normalization methods
# 1 - logp1
# 2 - CLR
# 3- Pearson correlation

# We will store normalization parameters (i.e. total counts) in the corresponding slot
# so when subsetting count matrix to the list of genes of interest (i.e. TFs) we will normalize them
# properly

# Parameters
# get total counts
get.SumCounts <- function(mat,
                          margin = 2L){
   if(!margin %in% c(1,2)){
      stop("Margin must be either 1 (cells as rows) or 2 (cells as column)")
   }
   
   if(margin == 2) {
   r <- Matrix::colSums(mat) 
   } else if (margin == 1) {
      r <- Matrix::rowSums(mat)
   }
   return (r)
}

GeomMean <- function(x){ 
   # compute geometric means
   #https://github.com/satijalab/seurat/blob/1549dcb3075eaeac01c925c4b4bb73c73450fc50/R/preprocessing.R#L4418
   exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) /
          length(x = x)
   )
}


get.GeomMean <- function(mat,
                         margin = 2L){
   if (margin == 1L) {
      # Geometric mean for rows
      gmeans <- sapply(1:nrow(mat), function(row_idx) {
         row_start <- mat@p[row_idx] + 1
         row_end <- mat@p[row_idx + 1]
         row_values <- mat@x[row_start:row_end]
         GeomMean(row_values)
      })
      names(gmeans) <- rownames(mat)  # Assign row names
   } else if (margin == 2L) {
      # Geometric mean for columns
      gmeans <- sapply(1:ncol(mat), function(col_idx) {
         col_start <- mat@p[col_idx] + 1
         col_end <- mat@p[col_idx + 1]
         col_values <- mat@x[col_start:col_end]
         GeomMean(col_values)
      })
      names(gmeans) <- colnames(mat)  # Assign column names
   } else {
      stop("Invalid margin. Use 1 for rows or 2 for columns.")
   }
   return(gmeans)   
}

# size factors for pearson residuals
# https://github.com/const-ae/transformGamPoi/blob/10bcccd4fc02e2659c2803ec97ccfbc6215e9600/R/helpers.R#L116
# size_factors <- transformGamPoi:::.handle_size_factors(size_factors = TRUE,
#                                                        mat,
#                                                        verbose = FALSE)

# join function to get params for normalization
get.Normalization_params <- function(mat,
                                     method = c("Log1p", "CLR", "pearson"),
                                     margin = 2L,
                                     size_factors = TRUE){
   # Run the requested methods
   norm_params <- switch(method[1],
                         "Log1p" = get.SumCounts(mat, margin),
                         "CLR" = get.GeomMean(mat, margin),
                         "pearson" = transformGamPoi:::.handle_size_factors(size_factors, Y = mat),
                         stop(method, " is not a supported normalization method.")
   )
   
   return(norm_params)
}


########################################################################################################
# Normalization 

Log_Normalize <- function(mat,
                          scale.factor = 1e4,
                          total.counts = NULL,
                          margin = 2L){
   
   if (is.null(total.counts)) {
      stop("`total.counts` must be provided")
   }
   
   cell.ids <- if (margin == 2L) colnames(mat) else rownames(mat)
   if (is.null(cell.ids) || !all(cell.ids %in% names(total.counts))) {
      stop("Missing normalization parameters: total counts per cell")
   }
   
   # Ensure counts are aligned
   total.counts <- total.counts[cell.ids]
   
   # Normalize matrix sparsely
   xnorm <- mat
   xnorm@x <- xnorm@x / rep(total.counts, diff(xnorm@p))  # Sparse division
   xnorm@x <- xnorm@x * scale.factor                     # Sparse scaling
   xnorm@x <- log1p(xnorm@x)  
   
   return(xnorm)
}



clr_Normalize <- function(mat,
                          GeomMeans = NULL,
                          margin = 2L){
   if (is.null(GeomMeans)) {
      stop("`GeomMeans` must be provided")
   }
   
   cell.ids <- if (margin == 2L) colnames(mat) else rownames(mat)
   if (is.null(cell.ids) || !all(cell.ids %in% names(GeomMeans))) {
      stop("Missing normalization parameters: Geometric means")
   }
   # Ensure geometric means are aligned
   GeomMeans <- GeomMeans[cell.ids]
   
   # CLR normalization sparingly
   xnorm <- mat
   xnorm@x <- xnorm@x / rep(GeomMeans, diff(xnorm@p))  # Sparse division
   xnorm@x <- log1p(xnorm@x)    
   
   return(xnorm)
}


# modified from https://github.com/const-ae/transformGamPoi/blob/10bcccd4fc02e2659c2803ec97ccfbc6215e9600/R/residual_transform.R#L101
residual_transform <- function(data,
                               residual_type = "pearson", #c("randomized_quantile", "pearson", "analytic_pearson"),
                               clipping = FALSE,
                               overdispersion = 0.05,
                               size_factors = NULL,
                               offset_model = TRUE,
                               overdispersion_shrinkage = TRUE,
                               ridge_penalty = 2,
                               on_disk = NULL,
                               return_fit = FALSE,
                               verbose = FALSE,
                               ...){
   
   if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
      message("Installing missing packages for pearson residuals normalization: glmGamPoi")
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
         message("Installing BiocManager...\n")
         install.packages("BiocManager")
      }
      BiocManager::install("glmGamPoi")
   }
   
   if (!requireNamespace("transformGamPoi", quietly = TRUE)) {
      stop("The 'transformGamPoi' package is required for pearson residuals normalization but not installed. Install it with: 
         devtools::install_github('const-ae/transformGamPoi')")
   }

   # tuned intput of size factors
   cell.ids <- colnames(data)
   if(!all(cell.ids %in% names(size_factors))){
      stop("Missing normalization parameters: Size factors")
   }
   
   # Get total counts per cell
   size_factors <- size_factors[cell.ids]
   
   
   # Allow any valid argument from glmGamPoi::residual.glmGamPoi()
   residual_type <- match.arg(residual_type[1], c("deviance", "pearson", "randomized_quantile",
                                                  "working", "response", "quantile", "analytic_pearson"))
   
   if(residual_type == "analytic_pearson"){
      return(transformGamPoi:::analytic_pearson_residual_transform(data = data,
                                                                   clipping = clipping,
                                                                   overdispersion = overdispersion,
                                                                   size_factors = size_factors,
                                                                   on_disk = on_disk,
                                                                   return_fit = return_fit,
                                                                   verbose = verbose))
   }
   
   if(inherits(data, "glmGamPoi")){
      fit <- data
   }else if(offset_model){
      counts <- transformGamPoi:::.handle_data_parameter(data, on_disk, allow_sparse = FALSE )
      
      fit <- glmGamPoi::glm_gp(counts, design = ~ 1,
                               size_factors = size_factors,
                               overdispersion = overdispersion,
                               overdispersion_shrinkage = overdispersion_shrinkage,
                               verbose = verbose,
                               ...)
   }else{
      counts <- transformGamPoi:::.handle_data_parameter(data, on_disk, allow_sparse = FALSE  )
      
      log_sf <- log(size_factors)
      attr(ridge_penalty, "target") <- c(0, 1)
      
      fit <- glmGamPoi::glm_gp(counts, design = ~ log_sf + 1,
                               size_factors = 1,
                               overdispersion = overdispersion,
                               overdispersion_shrinkage = overdispersion_shrinkage,
                               ridge_penalty = ridge_penalty,
                               verbose = verbose, ...)
   }
   
   if(overdispersion_shrinkage){
      # Use the dispersion trend when calculating the residuals
      fit$overdispersion_shrinkage_list$original_overdispersions <- fit$overdispersions
      fit$overdispersions <- fit$overdispersion_shrinkage_list$dispersion_trend
   }
   
   
   if(verbose){message("Calculate ", residual_type, " residuals")}
   
   resid <- stats::residuals(fit, type = residual_type)
   resid <- transformGamPoi:::clip_residuals(resid, clipping)
   resid <- transformGamPoi:::.convert_to_output(resid, data)
   
   if(! return_fit){
      resid
   }else{
      list(Residuals = resid, fit = fit)
   }
}



Normalize_data <- function(mat,
                           method = c("Log1p", "CLR", "pearson"),
                           margin = 2L,
                           scale.factor = 1e4,
                           size_factors = TRUE,
                           norm.params = NULL,
                           residual_type = "pearson",
                           ...){
   
   method <- method[1]
   
   if(method == "pearson"){
      mat <- as.matrix(mat)
   }
   
   if(is.null(norm.params)){
   norm.params <- get.Normalization_params(mat,
                                           method = method,
                                           margin = margin,
                                           size_factors = size_factors)
   }
   
   norm.mat <- switch(method,
                         "Log1p" = Log_Normalize(mat,
                                                 scale.factor = scale.factor,
                                                 total.counts = norm.params,
                                                 margin = margin),
                         "CLR" = clr_Normalize(mat,
                                               GeomMeans = norm.params,
                                               margin = margin),
                         "pearson" = residual_transform(data = mat,
                                                        size_factors = norm.params,
                                                        residual_type = residual_type,
                                                        ...),
                         stop(method, " is not a supported normalization method. Please use either Log1p, CLR, or pearson")
   )
   
}
