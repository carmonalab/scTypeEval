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
                         gene.lists = gene.lists,
                         black.list = black.list,
                         version = version)
   
   return(scTypeEval_obj)
}


set.activeIdent <- function(scTypeEval,
                            ident = NULL){
   if(is.null(ident)){
      stop("Specificy a cell type annotation in the provided metadata")
   }
   
   scTypeEval@active.ident <- ident
   
   return(scTypeEval)
}



add.HVG <- function(scTypeEval,
                    normalization.method = c("Log1p", "CLR", "pearson"),
                    ngenes = 500,
                    ...){
   # normalize matrix
   norm.mat <- Normalize_data(mat = scTypeEval@counts,
                              method = normalization.method[1],
                              ...)
   
   # get highly variable genes
   hgv <- get.HVG(norm.mat,
                  ngenes = ngenes)
   
   scTypeEval@gene.list[["HVG"]] <- hgv
   
   return(scTypeEval)
}


add.GeneMarkers <- function(scTypeEval,
                            method = c("scran.findMarkers", "gpsFISH"),
                            ngenes = 500){
   
   return(scTypeEval)
   
}


add.GeneList <- function(scTypeEval,
                         gene.list = NULL){
   
   # check if it is a list and if it is named
   # Check if gene.list is provided
   if (is.null(gene.list) || !is.list(gene.list)) {
      stop("gene.list cannot be NULL. Please provide a valid list.")
   }
   
   # If gene.list is unnamed, assign names
   if (any(is.null(names(gene.list)))) {
      names(gene.list) <- paste0("gene.list", seq_along(gene.list))
      warning("All or some names of the list is NULL, renaming list.")
   }
   
   scTypeEval@gene.lists <- c(scTypeEval@gene.lists, gene.list)
   
   return(scTypeEval)
}


# Wrapper for creating consistency assay

get.Consistency <- function(scTypeEval,
                            ident = NULL,
                            sample = NULL,
                            normalization.method = c("Log1p", "CLR", "pearson"),
                            gene.list = "all",
                            distance.method = "euclidean",
                            IntVal.metric = c("silhouette", "modularity", "ward",
                                              "inertia", "Xie-Beni", "S_Dbw", "I"),
                            data.type = c("sc", "pseudobulk", "pseudobulk_1vsall"),
                            ncores = 1,
                            bparam = NULL,
                            progressbar = TRUE,
                            ...
                            
){
   

   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   if(is.null(ident)){
      ident <- scTypeEval@active.ident
   }
   
   normalization.method <- normalization.method[1]
   data.type <- data.type[1]

   
   # get normalization parameters
   norm.params <- get.Normalization_params(mat = mat,
                                           method = normalization.method)

   # loop over each gene.list
   consist.list <- BiocParallel::bplapply(names(gene.list),
                                          BPPARAM = param,
                                          function(sc){
                                             # get data.type
                                             mat <- switch(data.type,
                                                           "sc" = scTypeEval@counts,
                                                           "pseudobulk" = get_pseudobulk())
                                             # subset gene list
                                             red.mat <- mat[gene.list[[sc]],]
                                             # normalize data
                                             norm.mat <- Normalize_data(red.mat,
                                                                        method = normalization.method,
                                                                        norm.params = norm.params)
                                             
                                             consist <- calculate_IntVal_metric(mat = red.mat,
                                                                                norm.mat = norm.mat,
                                                                                metrics = IntVal.metric,
                                                                                dist.method = distance.method,
                                                                                ident = ident,
                                                                                ...)
                                                
                                          })
   
   
   
}



