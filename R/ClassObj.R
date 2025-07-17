# S4 class object
methods::setClass("scTypeEval",
                  slots = c(
                     counts = "dgCMatrix", # raw counts
                     metadata = "data.frame", # data.frame with metadata
                     distances = "list",
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
   
   # Set default for consistency if not provided
   if (is.null(args$distances)) {
      args$distances <- list()
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
                             consistency = args$consistency, 
                             reductions = args$reductions)
   
   validObject(.Object)  # Validate the object
   .Object
})

methods::setClass("DissimilarityAssay",
                  slots = c(
                     measure = "ANY",
                     consistency.metric = "character",
                     distance.method = "ANY",
                     gene.list = "character",
                     black.list = "character",
                     ident = "character",
                     data.type = "character",
                     sample = "ANY"
                  )
)

# define consistency assay object
methods::setClass("ConsistencyAssay",
                  slots = c(
                     measure = "ANY",
                     consistency.metric = "character",
                     distance.method = "ANY",
                     gene.list = "character",
                     black.list = "character",
                     ident = "character",
                     data.type = "character",
                     sample = "ANY"
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
                     ident = "factor",
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
