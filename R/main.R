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
   scTypeEval_obj <- methods::new("scTypeEval",
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
   if(!ident %in% names(scTypeEval@metadata)){
      stop("Please provide a ident, i.e. a cell type or annotation
           to group cells included in metadata")
   }
   
   scTypeEval@active.ident <- ident
   
   return(scTypeEval)
}



add.HVG <- function(scTypeEval,
                    normalization.method = c("Log1p", "CLR", "pearson"),
                    ngenes = 500,
                    black.list = NULL,
                    ...){
   # normalize matrix
   mat <- scTypeEval@counts
   if(is.null(black.list)){
      black.list <- scTypeEval@black.list
   }
   # remove blacked listed genes
   mat <- mat[!rownames(mat) %in% black.list,]
   
   norm.mat <- Normalize_data(mat = mat,
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
                            min.samples = 5,
                            min.cells = 10,
                            black.list = NULL,
                            ncores = 1,
                            bparam = NULL,
                            progressbar = TRUE,
                            verbose = TRUE,
                            ...
                            
){
   if(is.null(ident)){
      ident <- scTypeEval@active.ident
   }
   
   if(!ident %in% names(scTypeEval@metadata)){
      stop("Please provide a ident, i.e. a cell type or annotation to group cells included in metadata")
   }
   
   # retrieve ident and convert to factor
   ident.name <- ident
   ident <- scTypeEval@metadata[[ident]]
   ident <- factor(ident)
   
   if(!is.null(sample)){
      if(!sample %in% names(scTypeEval@metadata)){
         stop("`sample` parameter not found in metadata colnames.")
      }
      # retrieve sample and convert to factor
      sample.name <- sample
      sample <- scTypeEval@metadata[[sample]]
      sample <- factor(sample)
      
   } else {
      if(data.type != "sc"){
         stop("For pseudobulk provide a dataset with multiple samples,
              and specifiy their respective column metadata in `sample` parameter")
      } else {
         if(verbose){message("Using dataset as a unique sample, computing consistency across cells.\n")}
      }
   }
   
   distance_methods <- c(
      "euclidean",
      "maximum",
      "manhattan",
      "canberra",
      "binary",
      "minkowski",
      "Jaccard",
      "Weighted_Jaccard",
      "gower",
      "bray-curtis",
      "cosine",
      "pearson"
   )
   if(distance.method %in% distance_methods){
      stop(distance, " distance method not supported. Pick up one of: ", 
           paste(distance_methods, collapse = ", "))
   }
   
   IntVal_metric = c("silhouette", "modularity", "ward",
                     "inertia", "Xie-Beni", "S_Dbw", "I")
   
   if(all(IntVal.metric %in% IntVal_metric)){
      stop(IntVal.metric, "at least one internal validation metrics(s) method not supported.
           Pick up one, some or all out of: ", 
           paste(IntVal_metric, collapse = ", "))
   }
   
   data_type = c("sc", "pseudobulk", "pseudobulk_1vsall")
   
   if(data.type %in% data_type){
      stop(data.type, " data type conversion method not supported. Pick up one of: ", 
           paste(data_type, collapse = ", "))
   }
   
   
   param1 <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   if(data_type == "pseudobulk_1vsall")
   

   
   normalization.method <- normalization.method[1]
   data.type <- data.type[1]
   
   if(is.null(black.list)){
   black.list <- scTypeEval@black.list
   }
   

   # loop over each gene.list
   consist.list <- BiocParallel::bplapply(names(gene.list),
                                          BPPARAM = param,
                                          function(sc){
                                             
                                             # get matrix
                                             mat <- scTypeEval@counts[gene.list[[sc]]]
                                             
                                             # remove black list genes
                                             mat <- mat[!rownames(mat) %in% black.list,]
                                             
                                             con.list <- consistency.helper(mat,
                                                                            ident = ident,
                                                                            sample = sample,
                                                                            normalization.method = normalization.method,
                                                                            distance.method = distance.method,
                                                                            IntVal.metric = IntVal.metric,
                                                                            data.type = data.type)
                                             
                                             # accommodte to ConsistencyAssay
                                             
                                            CA <- lapply(names(con.list),
                                                         function(cc){
                                                            methods::new("ConsistencyAssay",
                                                                         consistency.metric = con.list[[cc]],
                                                                         dist.method = distance.method,
                                                                         gene.list = sc,
                                                                         black.list = black.list,
                                                                         ident = ident.name,
                                                                         data.type = data.type,
                                                                         sample = sample.name)
                                                         })
                                            
                                            return(CA)
                                          })
   
   names(consist.list) <- names(gene.list)

   scTypeEval@consistency <- consist.list
   return(scTypeEval)

}


# get.ConsistencyData <- function(scTypeEval,
#                                 consisteny.metric = NULL,
#                                 dist.method = ){
#    
# }



