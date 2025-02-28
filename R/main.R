# Create scTypeEval object

create.scTypeEval <- function(matrix, 
                              metadata = NULL, 
                              gene.lists = list(), 
                              black.list = character(), 
                              active.ident = NULL,
                              version = "") {
   
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
   
   if(ncol(counts) != nrow(metadata)){
      stop("Different number of columns in counts and number of rows in metadata.")
   }
   
   # Create the scTypeEval object
   scTypeEval_obj <- methods::new("scTypeEval",
                                  counts = counts,
                                  metadata = metadata,
                                  gene.lists = gene.lists,
                                  black.list = black.list,
                                  active.ident = active.ident,
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
                    var.method = c("scran", "basic"),
                    sample = NULL,
                    ngenes = 500,
                    black.list = NULL,
                    ncores = 1,
                    bparam = NULL,
                    progressbar = TRUE,
                    ...){
   if(!is.null(sample)){
      if(!sample %in% names(scTypeEval@metadata)){
         stop("`sample` parameter not found in metadata colnames.")
      }
      # retrieve sample and convert to factor
      sample.name <- sample
      sample <- scTypeEval@metadata[[sample]]
      sample <- purge_label(sample)
      sample <- factor(sample)
   }
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   # normalize matrix
   mat <- scTypeEval@counts
   if(is.null(black.list)){
      black.list <- scTypeEval@black.list
   }
   
   
   norm.mat <- Normalize_data(mat = mat,
                              method = normalization.method[1],
                              ...)
   
   # remove blacked listed genes
   norm.mat <- norm.mat[!rownames(norm.mat) %in% black.list,]
   
   var.method <- var.method[1]
   
   # get highly variable genes
   hgv <- switch(var.method,
                 "basic" = get.HVG(norm.mat,
                                    ngenes = ngenes,
                                    sample = sample,
                                    bparam = param),
                 "scran" = get.GeneVar(norm.mat = norm.mat,
                                       sample = sample,
                                       ngenes = ngenes,
                                       bparam = param),
                 stop(var.method, " not supported for getting variable genes.")
                 )

   
   scTypeEval@gene.lists[["HVG"]] <- hgv
   
   return(scTypeEval)
}


add.GeneMarkers <- function(scTypeEval,
                            ident = NULL,
                            sample = NULL,
                            method = c("scran.findMarkers", "gpsFISH"),
                            ngenes.total = 500,
                            ngenes.celltype = 50,
                            ncores = 1,
                            bparam = NULL,
                            progressbar = TRUE,
                            ...){
   
   method <- method[1]
   if(!method %in% c("scran.findMarkers", "gpsFISH")){
      stop("Current gene markes definitions are either `scran.findMarkers` and `gpsFISH`")
   }
   
   if(is.null(ident)){
      ident <- scTypeEval@active.ident
   }
   
   if(!ident %in% names(scTypeEval@metadata)){
      stop("Please provide a ident, i.e. a cell type or annotation to group cells included in metadata")
   }
   
   # retrieve ident and convert to factor
   ident.name <- ident
   ident <- scTypeEval@metadata[[ident]]
   ident <- purge_label(ident)
   ident <- factor(ident)
   
   if(!is.null(sample)){
      if(!sample %in% names(scTypeEval@metadata)){
         stop("`sample` parameter not found in metadata colnames.")
      }
      # retrieve sample and convert to factor
      sample <- scTypeEval@metadata[[sample]]
      sample <- purge_label(sample)
      sample <- factor(sample)
   }
   
   if(is.null(black.list)){
      black.list <- scTypeEval@black.list
   }
   
   mat <- scTypeEval@counts
   
   markers <- switch(method,
                     "scran.findMarkers" = get.DEG(mat = mat,
                                                   ident = ident,
                                                   block = sample,
                                                   ngenes.celltype = ngenes.celltype,
                                                   ncores = ncores,
                                                   bparam = bparam,
                                                   progressbar = progressbar,
                                                   black.list = black.list,
                                                   ...),
                     "gpsFISH" = get.gpsFISHMarkers(sc_count = mat,
                                                    ident = ident,
                                                    black.list = black.list
                                                    )
                     )
   
   scTypeEval@gene.lists[[method]] <- markers
   
   
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

Run.Consistency <- function(scTypeEval,
                            ident = NULL,
                            sample = NULL,
                            normalization.method = c("Log1p", "CLR", "pearson"),
                            gene.list = NULL,
                            pca = FALSE,
                            ndim = 30,
                            distance.method = "euclidean",
                            IntVal.metric = c("silhouette", "NeighborhoodPurity",
                                              "ward.PropMatch", "Leiden.PropMatch",
                                              "modularity"),
                            data.type = c("sc", "pseudobulk",
                                          "pseudobulk_1vsall"),
                            min.samples = 5,
                            min.cells = 10,
                            KNNGraph_k = 5,
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
   
   data.type = data.type[1]
   
   # retrieve ident and convert to factor
   ident.name <- ident
   ident <- scTypeEval@metadata[[ident]]
   ident <- purge_label(ident)
   ident <- factor(ident)
   
   if(!is.null(sample)){
      if(!sample %in% names(scTypeEval@metadata)){
         stop("`sample` parameter not found in metadata colnames.")
      }
      # retrieve sample and convert to factor
      sample.name <- sample
      sample <- scTypeEval@metadata[[sample]]
      sample <- purge_label(sample)
      sample <- factor(sample)
      
      if(min.samples < 2){
         stop("For intersample comparison a minimum threshold of cell population
              in at least 2 samples is required, but more is recommended.")
      }
      
   } else {
      
      if(data.type != "sc"){
         stop("For pseudobulk provide a dataset with multiple samples,
              and specifiy their respective column metadata in `sample` parameter")
      } else {
         if(verbose){message("Using dataset as a unique sample, computing consistency across cells.\n")}
      }
      sample.name <- NULL
   }
   

   if(!distance.method %in% distance_methods){
      stop(distance.method, " distance method not supported. Pick up one of: ", 
           paste(distance_methods, collapse = ", "))
   }
   # only run EMD on pseudobulk data.type
   if(distance.method == "EMD" && data.type == "sc"){
      warning("Earth mover distance (EMD) for single-cell (sc) data.type is
              highly computially expensive... not recommended, it will take long")
   }
   
   
   if(!all(IntVal.metric %in% IntVal_metric)){
      stop("at least one internal validation metrics(s) method not supported.
           Pick up one, some or all out of: ", 
           paste(IntVal_metric, collapse = ", "))
   }
   
   if(!data.type %in% data_type){
      stop(data.type, " data type conversion method not supported. Pick up one of: ", 
           paste(data_type, collapse = ", "))
   }
   
   # set normalization method
   normalization.method <- normalization.method[1]
   
   # set gene lists
   if(is.null(gene.list)){
      gene.list <- scTypeEval@gene.lists
      if(length(gene.list) == 0){
         stop("No gene list found to run consistency metrics.\n
              Add custom gene list or compute highly variable genes with add.HVG()\n")
      }
   } else {
      if(!all(names(gene.list) %in% names(scTypeEval@gene.lists))){
         stop("Some gene list names not included in scTypeEval object")
      }
      gene.list <- scTypeEval@gene.lists[gene.list]
   }
   
   if(is.null(black.list)){
      black.list <- scTypeEval@black.list
   }
   
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   if(ncores == 1 && is.null(bparam)){
      param1 <- param
      param2 <- param
   } else {
      if(!is.null(sample) && data.type %in% c("sc", "pseudobulk_1vsall")){
         if(length(unique(sample))>length(gene.list)){
            param1 <- BiocParallel::SerialParam()
            param2 <- param
         } else {
            param1 <- param
            param2 <- BiocParallel::SerialParam()
         }
      } else {
         param1 <- param
         param2 <- BiocParallel::SerialParam()
      }
   }
   
  
   # loop over each gene.list
   consist.list <- BiocParallel::bplapply(names(gene.list),
                                          BPPARAM = param1,
                                          function(t){
                                             
                                             # get matrix
                                             keep <- rownames(scTypeEval@counts) %in% gene.list[[t]]
                                             mat <- scTypeEval@counts[keep,]

                                             # remove black list genes
                                             mat <- mat[!rownames(mat) %in% black.list,]
                                             
                                             con.list <- consistency.helper(mat,
                                                                            ident = ident,
                                                                            sample = sample,
                                                                            normalization.method = normalization.method,
                                                                            distance.method = distance.method,
                                                                            IntVal.metric = IntVal.metric,
                                                                            data.type = data.type,
                                                                            pca = pca,
                                                                            ndim = ndim,
                                                                            bparam = param2,
                                                                            min.samples = min.samples,
                                                                            min.cells = min.cells,
                                                                            KNNGraph_k = KNNGraph_k,
                                                                            verbose = verbose)
                                             
                                             # accommodte to ConsistencyAssay
                                             if(pca){
                                                t <- paste(t, "PCA", sep = ".")
                                             }
                                             
                                             CA <- lapply(seq_along(con.list),
                                                          function(cc){
                                                             lapply(names(con.list[[cc]]),
                                                                    function(y){
                                                                       if(!grepl(dist.need, y)){distance.method <- NULL}
                                                                       
                                                                       methods::new("ConsistencyAssay",
                                                                                    measure = con.list[[cc]][[y]],
                                                                                    consistency.metric = y,
                                                                                    distance.method = distance.method,
                                                                                    gene.list = t,
                                                                                    black.list = black.list,
                                                                                    ident = ident.name,
                                                                                    data.type = data.type,
                                                                                    sample = names(con.list)[[cc]])
                                                                       
                                                                    })
                                                             
                                                          }) %>% unlist()
                                            
                                    
                                             
                                             # name consistency assays
                                             names(CA) <- lapply(seq_along(CA),
                                                                 function(ca){
                                                                    v <- na.omit(c(CA[[ca]]@sample,
                                                                                   CA[[ca]]@data.type,
                                                                                   CA[[ca]]@distance.method,
                                                                                   CA[[ca]]@consistency.metric))
                                                                    paste(v,
                                                                          collapse = "_")
                                                                 })
                                            
                                            return(CA)
                                          })
   if(pca){
      names(consist.list) <- paste(names(gene.list), "PCA", sep = ".")
   } else {
      names(consist.list) <- names(gene.list)
   }
   
   # add to scTypeEval object
   for(n in names(consist.list)){
      for(m in names(consist.list[[n]]))
         scTypeEval@consistency[[n]][[m]] <- consist.list[[n]][[m]]
   }
   
   return(scTypeEval)

}


Run.BestHit <- function(scTypeEval,
                        ident = NULL,
                        ident_GroundTruth = NULL,
                        method = c("Mutual.Score", "Mutual.Match"),
                        sample = NULL,
                        data.type = "sc",
                        gene.list = NULL,
                        black.list = NULL,
                        min.cells = 10,
                        min.samples = 5,
                        ncores = 1,
                        bparam = NULL,
                        progressbar = TRUE,
                        verbose = TRUE,
                        ...
                            
){
   
   if(!data.type %in% data_type){
      stop(data.type, " data type conversion method not supported. Pick up one of: ", 
           paste(data_type[1:2], collapse = ", "))
   } else if(data.type == "pseudobulk_1vsall"){
      stop("pseudobulk_1vsall not supported for mutual BestHit.") 
   }
   
   if(is.null(ident)){
      ident <- scTypeEval@active.ident
   }
   
   if(!ident %in% names(scTypeEval@metadata)){
      stop("Please provide a ident, i.e. a cell type or annotation to group cells included in metadata")
   }

   # retrieve ident and convert to factor
   ident.name <- ident
   ident <- scTypeEval@metadata[[ident]]
   ident <- purge_label(ident)
   ident <- factor(ident)
   
   if(!is.null(ident_GroundTruth)){
      ident_GroundTruth <- scTypeEval@metadata[[ident_GroundTruth]]
      ident_GroundTruth <- purge_label(ident_GroundTruth)
      ident_GroundTruth <- factor(ident_GroundTruth)
   }
   
   if(!is.null(sample)){
      if(!sample %in% names(scTypeEval@metadata)){
         stop("`sample` parameter not found in metadata colnames.")
      }
      # retrieve sample and convert to factor
      sample.name <- sample
      sample <- scTypeEval@metadata[[sample]]
      sample <- purge_label(sample)
      sample <- factor(sample)
   } else {
      stop("BestHit consistency requires multiple samples and specify it in `sample` parameter.")
   }

   # set gene lists
   if(is.null(gene.list)){
      gene.list <- scTypeEval@gene.lists
      if(length(gene.list) == 0){
         stop("No gene list found to run consistency metrics.\n
              Add custom gene list or compute highly variable genes with add.HVG()\n")
      }
   } else {
      if(!all(names(gene.list) %in% names(scTypeEval@gene.lists))){
         stop("Some gene list names not included in scTypeEval object")
      }
      gene.list <- scTypeEval@gene.lists[gene.list]
   }
   
   if(is.null(black.list)){
      black.list <- scTypeEval@black.list
   }
   
   # set methods
   if(any(!method %in% mutual_method)){
      stop("Not supported consistency metrics, please provide both or either: ", mutual_method)
   }
   
   if(data.type == "sc" && "Mutual.Match" %in% method){
      warning("Mutual.Match consistency only supported for pseudobulk data.type, not running")
      method <- method[method != "Mutual.Match"]
   }
   
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)

   # loop over each gene.list
   consist.list <- lapply(names(gene.list),
                          function(t){
                             
                             # get matrix
                             keep <- rownames(scTypeEval@counts) %in% gene.list[[t]]
                             mat <- scTypeEval@counts[keep,]
                             
                             # remove black list genes
                             mat <- mat[!rownames(mat) %in% black.list,]
                             
                             con <- bestHit.SingleR(mat = mat,
                                                    ident = ident,
                                                    ident_GroundTruth = ident_GroundTruth,
                                                    sample = sample,
                                                    data.type = data.type,
                                                    method = method,
                                                    min.cells = min.cells,
                                                    min.samples = min.samples,
                                                    bparam = param)
                             
                             # accommodte to ConsistencyAssay
                             
                             CA <- lapply(names(con),
                                          function(cc){
                                             
                                             methods::new("ConsistencyAssay",
                                                          measure = con[[cc]],
                                                          consistency.metric = paste("BestHit", cc, sep = "-"),
                                                          distance.method = NA,
                                                          gene.list = t,
                                                          black.list = black.list,
                                                          ident = ident.name,
                                                          data.type = data.type,
                                                          sample = NA)
                                             
                                          })
                             
                             # name consistency assays
                             names(CA) <- lapply(seq_along(CA),
                                                 function(ca){
                                                    v <- na.omit(c(CA[[ca]]@data.type,
                                                                   CA[[ca]]@consistency.metric))
                                                    paste(v,
                                                          collapse = "_")
                                                 })
                             
                             
                             return(CA)
                          })
   
   names(consist.list) <- names(gene.list)
   
   # add to scTypeEval object
   for(n in names(consist.list)){
      for(m in names(consist.list[[n]]))
         scTypeEval@consistency[[n]][[m]] <- consist.list[[n]][[m]]
   }
   
   return(scTypeEval)
   
}

## wrapper of consistency metrics and BestHit

Run.scTypeEval <- function(scTypeEval,
                           ident = NULL,
                           ident_GroundTruth = NULL,
                           sample = NULL,
                           normalization.method = c("Log1p", "CLR", "pearson"),
                           gene.list = NULL,
                           pca = FALSE,
                           ndim = 30,
                           distance.method = "euclidean",
                           IntVal.metric = c("silhouette", "NeighborhoodPurity",
                                             "ward.PropMatch", "Leiden.PropMatch",
                                             "modularity"),
                           BH.method = c("Mutual.Score", "Mutual.Match"),
                           data.type = c("sc", "pseudobulk",
                                         "pseudobulk_1vsall"),
                           min.samples = 5,
                           min.cells = 10,
                           k.sc = 20,
                           k.psblk = 5,
                           black.list = NULL,
                           ncores = 1,
                           bparam = NULL,
                           progressbar = FALSE,
                           verbose = TRUE){

   
   df.res <- list()
   
   # only run EMD on pseudobulk data.type
   if(distance.method == "EMD" && "sc" %in% data.type){
      warning("Earth mover distance (EMD) for single-cell (sc) data.type is
              highly computially expensive... not running consistency for sc.\n")
      data.type <- data.type[data.type != "sc"]
   }
   
   if(length(data.type) == 0 | is.null(data.type)){
      stop("No data.type provided to run consistency.")
   }
   
   for(dt in data.type){
      
      if(dt == "sc"){
         spl = NULL
         k <- k.sc
      } else {
         spl = sample
         k <- k.psblk
      }
      
      if(verbose){message("\n------- Running Consistency for ", dt, " ",  format(Sys.time()), "\n")}
      
      sc.tmp <- Run.Consistency(scTypeEval,
                                ident = ident,
                                sample = spl,
                                normalization.method = normalization.method,
                                distance.method = distance.method,
                                IntVal.metric = IntVal.metric,
                                gene.list = gene.list,
                                data.type =  dt,
                                pca = pca,
                                ndim = ndim,
                                min.samples = min.samples,
                                min.cells = min.cells,
                                black.list = black.list,
                                ncores = ncores,
                                KNNGraph_k = k,
                                progressbar = progressbar,
                                verbose = FALSE)
      
      if(dt != "pseudobulk_1vsall" && any(BH.method %in% mutual_method)){
         
         if(verbose){message("------- Running BestHit for ", dt, " ",  format(Sys.time()), "\n")}
         sc.tmp <- Run.BestHit(sc.tmp,
                               data.type = dt,
                               ident = ident,
                               ident_GroundTruth = ident_GroundTruth,
                               method = BH.method,
                               gene.list = gene.list,
                               sample = sample,
                               min.samples = min.samples,
                               min.cells = min.cells,
                               black.list = black.list,
                               ncores = ncores,
                               progressbar = progressbar,
                               verbose = FALSE)
      }
      
      df <- get.ConsistencyData(sc.tmp)
      df.res[[dt]] <- df
   }
   
   wdf <- do.call(rbind, df.res)
   
   return(wdf)
}


get.ConsistencyData <- function(scTypeEval,
                                gene.list = NULL,
                                consistency.metric = NULL,
                                distance.method = NULL,
                                data.type = NULL
)
{
   
   # Ensure the consistency slot is not empty
   if (length(scTypeEval@consistency) == 0) {
      stop("The 'consistency' slot in the 'scTypeEval' object is empty.")
   }
   
   # Initialize an empty list to store filtered data
   filtered_data <- list()
   assays <- unlist(scTypeEval@consistency)
   
   # Loop through each ConsistencyAssay object in scTypeEval@consistency
   for (a in seq_along(assays)) {
      assay <- assays[[a]]
      # Check for class validity
      if (!inherits(assay, "ConsistencyAssay")) {
         stop(paste("Invalid object in consistency slot"))
      }
      
      # Apply filtering:
      # For `gene.list`, check if it is NULL or if the intersection with the assay's gene list is non-empty
      if (!is.null(gene.list) && length(intersect(gene.list, assay@gene.list)) == 0) next
      # For other parameters, check if they are NULL or match any value in the corresponding vector
      if (!is.null(consistency.metric) && !assay@consistency.metric %in% consistency.metric) next
      if (!is.null(distance.method) && !identical(assay@distance.method, distance.method)) next
      if (!is.null(data.type) && !assay@data.type %in% data.type) next
      
      cm <- assay@consistency.metric
      
      # scale if needed, convert to all metric to a -1 to 1 scale
      scaled_measure <- normalize_metric(value = assay@measure,
                                         metric = cm)
     
      
      # Extract relevant information into a named list
      filtered_data[[a]] <- data.frame(
         celltype = names(assay@measure),
         measure = assay@measure,
         scaled_measure = scaled_measure,
         consistency.metric = cm,
         distance.method = if (is.null(assay@distance.method)) NA else assay@distance.method,
         gene.list = paste(assay@gene.list, collapse = ", "),
         ident = assay@ident,
         data.type = assay@data.type,
         sample = if (is.null(assay@sample)) NA else assay@sample
      )
   }
   
   # Combine all data frames into one
   if (length(filtered_data) == 0) {
      stop("Filters provided yielded no result in current scTypeEval object") # Return an empty data frame if no results
   }
   
   result <- do.call(rbind, filtered_data)

   
   return(result)
   
   
}



# obtain PCA from a gene list
add.PCA <- function(scTypeEval,
                    ident = NULL,
                    sample = NULL,
                    normalization.method = c("Log1p", "CLR", "pearson"),
                    gene.list = NULL,
                    data.type = c("sc", "pseudobulk",
                                  "pseudobulk_1vsall"),
                    min.samples = 5,
                    min.cells = 10,
                    black.list = NULL,
                    ndim = 30,
                    ncores = 1,
                    bparam = NULL,
                    progressbar = TRUE,
                    verbose = TRUE){
   
   if(is.null(ident)){
      ident <- scTypeEval@active.ident
   }
   
   if(!ident %in% names(scTypeEval@metadata)){
      stop("Please provide a ident, i.e. a cell type or annotation to group cells included in metadata")
   }
   
   data.type = data.type[1]
   
   # retrieve ident and convert to factor
   ident.name <- ident
   ident <- scTypeEval@metadata[[ident]]
   ident <- purge_label(ident)
   ident <- factor(ident)
   
   if(!is.null(sample)){
      if(!sample %in% names(scTypeEval@metadata)){
         stop("`sample` parameter not found in metadata colnames.")
      }
      # retrieve sample and convert to factor
      sample.name <- sample
      sample <- scTypeEval@metadata[[sample]]
      sample <- purge_label(sample)
      sample <- factor(sample)
      
   } else {
      
      if(data.type != "sc"){
         stop("For pseudobulk provide a dataset with multiple samples,
              and specifiy their respective column metadata in `sample` parameter")
      } else {
         if(verbose){message("Using dataset as a unique sample, computing PCA across cells.\n")}
      }
      sample.name <- NULL
   }
   
   if(!data.type %in% data_type){
      stop(data.type, " data type conversion method not supported. Pick up one of: ", 
           paste(data_type, collapse = ", "))
   }
   
   # set normalization method
   normalization.method <- normalization.method[1]
   
   # set gene lists
   if(is.null(gene.list)){
      gene.list <- scTypeEval@gene.lists
      if(length(gene.list) == 0){
         stop("No gene list found to run PCA.\n
              Add custom gene list or compute highly variable genes with add.HVG()\n")
      }
   } else {
      if(!all(names(gene.list) %in% names(scTypeEval@gene.lists))){
         stop("Some gene list names not included in scTypeEval object")
      }
      gene.list <- scTypeEval@gene.lists[gene.list]
   }
   
   if(is.null(black.list)){
      black.list <- scTypeEval@black.list
   }
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   # loop over each gene.list
   pca.list <- BiocParallel::bplapply(names(gene.list),
                                          BPPARAM = param,
                                          function(t){
                                             
                                             # get matrix
                                             keep <- rownames(scTypeEval@counts) %in% gene.list[[t]]
                                             mat <- scTypeEval@counts[keep,]
                                             
                                             # remove black list genes
                                             mat <- mat[!rownames(mat) %in% black.list,]
                                             
                                             # compute PCA
                                             pc.list <- get.PCA(mat,
                                                                ident = ident,
                                                                sample = sample,
                                                                normalization.method = normalization.method,
                                                                data.type = data.type,
                                                                min.samples = min.samples,
                                                                min.cells = min.cells)
                                             
                                             # accommodte to DimReduc
                                             
                                             pc.list <- lapply(seq_along(pc.list),
                                                               function(cc){
                                                                  pc.list[[cc]]@gene.list <- t
                                                                  pc.list[[cc]]@black.list <- black.list
                                                                  pc.list[[cc]]@sample <- names(pc.list)[[cc]]
                                                                  return(pc.list[[cc]])
                                                               })
                                             
                                             
                                             
                                             # name DimRed assays
                                             names(pc.list) <- lapply(seq_along(pc.list),
                                                                 function(ca){
                                                                    v <- na.omit(c(pc.list[[ca]]@sample,
                                                                                   pc.list[[ca]]@data.type))
                                                                    paste(v,
                                                                          collapse = "_")
                                                                 })
                                             
                                             return(pc.list)
                                          })
   
   names(pca.list) <- names(gene.list)
   
   # add to scTypeEval object
   for(n in names(pca.list)){
      for(m in names(pca.list[[n]]))
         scTypeEval@reductions[[n]][[m]] <- pca.list[[n]][[m]]
   }
   
   return(scTypeEval)
   
}



plot.PCA <- function(scTypeEval,
                     gene.list = NULL,
                     data.type = NULL,
                     dims = c(1,2),
                     label = TRUE,
                     show.legend = F
                     ){
   
   # Ensure the reductions slot is not empty
   if (length(scTypeEval@reductions) == 0) {
      stop("The 'reductions' slot in the 'scTypeEval' object is empty.")
   }
   
   assays <- unlist(scTypeEval@reductions)
   
   # Loop through each ConsistencyAssay object in scTypeEval@consistency
   pls <- lapply(names(assays),
                 function(a) {
                    assay <- assays[[a]]
                    # Check for class validity
                    if (!inherits(assay, "DimRed")) {
                       stop(paste("Invalid object in reductions slot"))
                    }
                    
                    # Apply filtering:
                    # For `gene.list`, check if it is NULL or if the intersection with the assay's gene list is non-empty
                    if (!is.null(gene.list) && length(intersect(gene.list, assay@gene.list)) == 0) next
                    if (!is.null(data.type) && !assay@data.type %in% data.type) next
                    
                    
                    df <- assay@embeddings[,dims] |>
                       as.data.frame() |>
                       dplyr::mutate(ident = assay@ident)
                    labs <- paste0("PC", dims)
                    
                    pl <- helper.plot.PCA(df,
                                          show.legend = show.legend,
                                          label = label) +
                       ggplot2::labs(x = labs[1],
                                     y = labs[2],
                                     title = a)
                    
                    return(pl)
                 })
   
   names(pls) <- names(assays)
   
   return(pls)
   
}
