




# join function to get params for normalization
add.Normalization_params <- function(scTypeEval,
                                     method = c("Log1p", "CLR", "pearson"),
                                     margin = 2L,
                                     size_factors = TRUE,
                                     return_value = F){
   
   mat <- scTypeEval@counts
   
   # Run the requested methods
   norm_params <- switch(method[1],
                         "Log1p" = get.SumCounts(mat, margin),
                         "CLR" = get.GeoMean(dist, ident, KNNGraph_k),
                         "pearson" = transformGamPoi:::.handle_size_factors(size_factors, mat),
                         stop(metric, " is not a supported normalization method.")
   )
   
   if(return_value){
      return(norm_params)
   } else {
      scTypeEval@norm.param[[method]] <- norm_params
      return(scTypeEval)
   }
   
}
