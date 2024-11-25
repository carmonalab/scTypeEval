# we will offer 3 types of normalization methods
# 1 - logp1
# 2 - CLR
# 3- Pearson correlation

# We will store normalization parameters (i.e. total counts) in the corresponding slot
# so when subsetting count matrix to the list of genes of interest (i.e. TFs) we will normalize them
# properly

# Parameters
# get total counts
get.SumCounts <- function(mat){
   colSums(mat)
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
   gmeans <- apply(mat, MARGIN = margin, GeomMean)
   return(gmeans)   
}

# size factors for pearson residuals
# https://github.com/const-ae/transformGamPoi/blob/10bcccd4fc02e2659c2803ec97ccfbc6215e9600/R/helpers.R#L116
size_factors <- transformGamPoi:::.handle_size_factors(size_factors = TRUE,
                                                       mat,
                                                       verbose = FALSE)



# Normalization

Log_Normalize <- function(mat,
                          scale.factor = 1e4,
                          total.counts = NULL,
                          margin = 2L){
   
   cell.ids <- colnames(mat)
   if(!all(cell.ids %in% names(total.counts))){
      stop("Missing normalization parameters: total counts per cell")
   }
   
   # Get total counts per cell
   total.counts <- total.counts[cell.ids]
   
   # Normalize the matrix by scaling the values using the precomputed total counts per cell
   xnorm <- sweep(mat, MARGIN = 2, total.counts, FUN = "/")  # Divide by total counts per cell
   # Scale by the given scale factor
   xnorm <- sweep(xnorm, MARGIN = 2, scale.factor, FUN = "*")
   # Apply log1p normalization efficiently
   xnorm <- log1p(xnorm)
   
   return(xnorm)
}



clr_Normalize <- function(mat,
                          GeoMeans = NULL){
   cell.ids <- colnames(mat)
   if(!all(cell.ids %in% names(GeoMeans))){
      stop("Missing normalization parameters: Geometric mean")
   }
   
   # Get total counts per cell
   GeomMeans <- GeoMeans[cell.ids]
   
   # Divide by geometric means
   xnorm <- sweep(mat, MARGIN = 2, GeoMeans, FUN = "/")  
   # Apply log1p normalization efficiently
   xnorm <- log1p(xnorm)
   
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
