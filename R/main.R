# Create scTypeEval object

create.scTypeEval <- function(matrix, 
                              metadata = NULL, 
                              gene.lists = list(), 
                              black.list = character(), 
                              active.ident = NULL) {
   
   # Check input type
   if (inherits(matrix, "Seurat")) {
      counts <- as(matrix@assays$RNA@counts, "dgCMatrix")
      metadata <- as.data.frame(matrix@meta.data)
   } else if (inherits(matrix, "SingleCellExperiment")) {
      counts <- as(assay(matrix, "counts"), "dgCMatrix")
      metadata <- as.data.frame(colData(matrix))
   } else if (inherits(matrix, "matrix") || inherits(matrix, "dgCMatrix")) {
      counts <- as(matrix, "dgCMatrix")
      if (is.null(metadata)) {
         stop("For matrix input, metadata dataframe must be provided.")
      }
      metadata <- as.data.frame(metadata)
   } else {
      stop("Input object must be a Seurat, SingleCellExperiment, or matrix-like object.")
   }
   
   # Create the scTypeEval object
   scTypeEval_obj <- new("scTypeEval",
                         counts = counts,
                         metadata = metadata,
                         norm.param = norm.param,
                         gene.lists = gene.lists,
                         black.list = black.list,
                         version = version)
   
   return(scTypeEval_obj)
}





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
