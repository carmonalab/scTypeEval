# S4 class object
methods::setClass("scTypeEval",
                  slots = c(
                     counts = "dgCMatrix", # raw counts
                     metadata = "data.frame", # data.frame with metadata
                     consistency = "list", # actual consistency results assays
                     gene.lists = "list", # list of either HGV, markers
                     black.list = "character", # list of genes in black list
                     active.ident = "factor", # default grouping variable
                     reductions = "list", # list of dim reductions
                     misc = "list", # miscellaneous
                     commands = "list", # commands performed
                     version = "character" # package version
                  )
)

# Define the initialize method to set default for reductions and consistency
methods::setMethod("initialize", "scTypeEval", function(.Object, ...) {
   args <- list(...)
   
   # Set default for consistency if not provided
   if (is.null(args$consistency)) {
      args$consistency <- list()
   }
   
   # Set default for reductions if not provided
   if (is.null(args$reductions)) {
      args$reductions <- list()
   }
   
   # Set default for metadata if not provided
   if (is.null(args$metadata)) {
      args$metadata <- data.frame()
   }
   
   # Pass the updated arguments to the default initialize method
   .Object <- callNextMethod(.Object, ..., 
                             consistency = args$consistency, 
                             reductions = args$reductions, 
                             metadata = args$metadata)
   
   validObject(.Object)  # Validate the object
   .Object
})

# define consistency assay object
methods::setClass("ConsistencyAssay",
                  slots = c(
                     metrics = "data.frame",
                     method = "character",
                     gene.list = "character",
                     black.list = "character",
                     ident = "factor"
                  )
)

methods::setClass("DimReduc",
                  slots = c(
                     cell.embeddings = 'matrix',
                     feature.loadings = 'matrix',
                     feature.loadings.projected = 'matrix',
                     gene.list = "character",
                     black.list = "character",
                     key = "character" # type of reduction, PCA, UMAP...
                  )
)
