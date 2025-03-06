hirearchy.helper <- function(mat,
                             ident,
                             normalization.method,
                             distance.method = "euclidean",
                             hirearchy.method = "ward.D2",
                             data.type,
                             sample = NULL,
                             pca = FALSE,
                             ndim = 30,
                             min.samples = 5,
                             min.cells = 10,
                             verbose = TRUE){
   
   mat <- get.matrix(mat,
                     data.type = data.type,
                     ident = ident,
                     sample = sample,
                     min.samples = min.samples,
                     min.cells = min.cells,
                     bparam = bparam)
   
   norm.mat <- Normalize_data(mat@matrix,
                              method = normalization.method)
   if(!pca){
      mat <- mat@matrix
   } else {
      pr <- custom_prcomp(norm.mat, ndim)
      mat <- NULL
      norm.mat <- t(pr$x)
   }
   
   dist <- get.distance(mat = mat,
                        norm.mat = norm.mat,
                        distance.method = distance.method)
   hclust_result <- stats::hclust(dist,
                                  method = hirearchy.method)
   
   return(hclust_result)
   
}
