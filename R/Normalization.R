# we will offer 3 types of normalization methods
# 1 - logp1
# 2 - CLR
# 3- Pearson correlation

# We will store normalization parameters (i.e. total counts) in the corresponding slot
# so when subsetting count matrix to the list of genes of interest (i.e. TFs) we will normalize them
# properly

# Parameters
# get total counts
get_sum_counts <- function(mat,
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

geom_mean <- function(x){ 
   # compute geometric means
   #https://github.com/satijalab/seurat/blob/1549dcb3075eaeac01c925c4b4bb73c73450fc50/R/preprocessing.R#L4418
   exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) /
          length(x = x)
   )
}


get_geom_mean <- function(mat,
                         margin = 2L){
   if (margin == 1L) {
      # Geometric mean for rows
      gmeans <- vapply(seq_len(nrow(mat)), function(row_idx) {
         row_start <- mat@p[row_idx] + 1
         row_end <- mat@p[row_idx + 1]
         row_values <- mat@x[row_start:row_end]
         geom_mean(row_values)
      }, numeric(1))
      names(gmeans) <- rownames(mat)
   } else if (margin == 2L) {
      # Geometric mean for columns
      gmeans <- vapply(seq_len(ncol(mat)), function(col_idx) {
         col_start <- mat@p[col_idx] + 1
         col_end <- mat@p[col_idx + 1]
         col_values <- mat@x[col_start:col_end]
         geom_mean(col_values)
      }, numeric(1))
      names(gmeans) <- colnames(mat)
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
get_normalization_params <- function(mat,
                                     method = c("Log1p", "CLR", "pearson"),
                                     margin = 2L,
                                     size_factors = TRUE){
   if(method[1] == "pearson"){
      # Check required packages only when pearson is requested
      check_residuals_pck()
   }
   # Run the requested methods
   norm_params <- switch(method[1],
                         "Log1p" = get_sum_counts(mat, margin),
                         "CLR" = get_geom_mean(mat, margin),
                         "pearson" = getFromNamespace(".handle_size_factors", "transformGamPoi")(size_factors, Y = mat),
                         stop(method, " is not a supported normalization method.")
   )
   
   # Assign names if missing
   if (is.null(names(norm_params))) {
      names(norm_params) <- if (margin == 2L) colnames(mat) else rownames(mat)
   }
   
   return(norm_params)
}


# Normalization 

log_normalize <- function(mat,
                          scale_factor = 1e4,
                          total_counts = NULL,
                          margin = 2L){
   
   if (is.null(total_counts)) {
      stop("`total_counts` must be provided")
   }
   
   cell.ids <- if (margin == 2L) colnames(mat) else rownames(mat)
   if (is.null(cell.ids) || !all(cell.ids %in% names(total_counts))) {
      stop("Missing normalization parameters: total counts per cell")
   }
   
   # Ensure counts are aligned
   total_counts <- total_counts[cell.ids]
   
   # Normalize matrix sparsely
   xnorm <- mat
   xnorm@x <- xnorm@x / rep(total_counts, diff(xnorm@p))  # Sparse division
   xnorm@x <- xnorm@x * scale_factor                     # Sparse scaling
   xnorm@x <- log1p(xnorm@x)  
   
   return(xnorm)
}



clr_normalize <- function(mat,
                          geom_means = NULL,
                          margin = 2L){
   if (is.null(geom_means)) {
      stop("`geom_means` must be provided")
   }
   
   cell.ids <- if (margin == 2L) colnames(mat) else rownames(mat)
   if (is.null(cell.ids) || !all(cell.ids %in% names(geom_means))) {
      stop("Missing normalization parameters: Geometric means")
   }
   # Ensure geometric means are aligned
   geom_means <- geom_means[cell.ids]
   
   # CLR normalization sparingly
   xnorm <- mat
   xnorm@x <- xnorm@x / rep(geom_means, diff(xnorm@p))  # Sparse division
   xnorm@x <- log1p(xnorm@x)    
   
   return(xnorm)
}

check_residuals_pck <- function(){
   # Only used when pearson normalization is requested; fail fast if missing
   required_packages <- c("transformGamPoi", "glmGamPoi")
   missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
   if (length(missing_packages) > 0) {
      stop("Missing required packages for pearson residuals normalization: ",
           paste(missing_packages, collapse = ", "),
           "\nPlease install them before using method = 'pearson'.", call. = FALSE)
   }
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
   # Check and install required packages if missing
   check_residuals_pck()

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
      tg_analytic <- getFromNamespace("analytic_pearson_residual_transform", "transformGamPoi")
      return(tg_analytic(data = data,
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
      tg_handle_data <- getFromNamespace(".handle_data_parameter", "transformGamPoi")
      counts <- tg_handle_data(data, on_disk, allow_sparse = FALSE )
      
      fit <- glmGamPoi::glm_gp(counts, design = ~ 1,
                               size_factors = size_factors,
                               overdispersion = overdispersion,
                               overdispersion_shrinkage = overdispersion_shrinkage,
                               verbose = verbose,
                               ...)
   }else{
      tg_handle_data <- getFromNamespace(".handle_data_parameter", "transformGamPoi")
      counts <- tg_handle_data(data, on_disk, allow_sparse = FALSE  )
      
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
   tg_clip <- getFromNamespace("clip_residuals", "transformGamPoi")
   tg_convert <- getFromNamespace(".convert_to_output", "transformGamPoi")
   resid <- tg_clip(resid, clipping)
   resid <- tg_convert(resid, data)
   
   if(! return_fit){
      resid
   }else{
      list(Residuals = resid, fit = fit)
   }
}



normalize_data <- function(mat,
                           method = c("Log1p", "CLR", "pearson"),
                           margin = 2L,
                           scale_factor = 1e4,
                           size_factors = TRUE,
                           norm_params = NULL,
                           residual_type = "pearson",
                           ...){
   
   method <- method[1]
   
   if(method == "pearson"){
      mat <- as.matrix(mat)
   }
   
   if(is.null(norm_params)){
   norm_params <- get_normalization_params(mat,
                                           method = method,
                                           margin = margin,
                                           size_factors = size_factors)
   }
   
   norm_mat <- switch(method,
                         "Log1p" = log_normalize(mat,
                                                 scale_factor = scale_factor,
                                                 total_counts = norm_params,
                                                 margin = margin),
                         "CLR" = clr_normalize(mat,
                                               geom_means = norm_params,
                                               margin = margin),
                         "pearson" = residual_transform(data = mat,
                                                        size_factors = norm_params,
                                                        residual_type = residual_type,
                                                        ...),
                         stop(method, " is not a supported normalization method. Please use either Log1p, CLR, or pearson")
   )
   
   # Convert to sparse matrix if needed (pearson returns dense matrix)
   if (method == "pearson" && !inherits(norm_mat, "dgCMatrix")) {
      norm_mat <- Matrix::Matrix(norm_mat, sparse = TRUE)
   }
   
   return(norm_mat)
}
