# S4 class object
methods::setClass("scTypeEval",
                  slots = c(
                     counts = "dgCMatrix", # raw counts
                     metadata = "data.frame", # data.frame with metadata
                     data = "list",
                     dissimilarity = "list",
                     consistency = "list", # actual consistency results assays
                     gene.lists = "list", # list of either HGV, markers
                     black.list = "character", # list of genes in black list
                     active.ident = "ANY", # default grouping variable
                     reductions = "list", # list of dim reductions
                     misc = "list", # miscellaneous
                     version = "character" # package version
                  )
)

# Define the initialize method to set default for reductions and consistency
methods::setMethod("initialize", "scTypeEval", function(.Object, ...) {
   args <- list(...)
   
   if (is.null(args$data)) {
      args$data <- list()
   }
   
   # Set default for consistency if not provided
   if (is.null(args$dissimilarity)) {
      args$dissimilarity <- list()
   }
   
   # Set default for consistency if not provided
   if (is.null(args$consistency)) {
      args$consistency <- list()
   }
   
   # Set default for reductions if not provided
   if (is.null(args$reductions)) {
      args$reductions <- list()
   }
   
   # Pass the updated arguments to the default initialize method
   .Object <- callNextMethod(.Object, ..., 
                             data = args$data,
                             dissimilarity = args$dissimilarity,
                             consistency = args$consistency, 
                             reductions = args$reductions)
   
   validObject(.Object)  # Validate the object
   .Object
})

methods::setClass("DataAssay",
                  slots = c(
                     data = 'matrix',
                     aggregation = "character",
                     group = "factor",
                     sample = "factor",
                     ident = "list"
                  )
)

methods::setClass("DissimilarityAssay",
                  slots = c(
                     dissimilarity = "dist",
                     method = "character",
                     aggregation = "character",
                     gene.list = "character",
                     black.list = "character",
                     ident = "list",
                     sample = "factor"
                  )
)


methods::setClass("DimRed",
                  slots = c(
                     embeddings = 'matrix',
                     feature.loadings = 'matrix',
                     gene.list = "character",
                     black.list = "character",
                     aggregation = "character",
                     group = "factor",
                     sample = "factor",
                     ident = "list",
                     key = "character" # type of reduction, PCA, UMAP...
                  )
)

methods::setClass("Mat_ident",
                  slots = c(
                     matrix = 'matrix',
                     groups = 'factor',
                     ident = 'factor',
                     sample = 'factor'
                  )
)
