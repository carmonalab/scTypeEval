# auxiliary function to process all dataets, without filtering for obtaining the PCA or 
# matrix with relevant genes on the gene list

process_clustering <- function(scTypeEval,
                               sample = "sample",
                               reduction = TRUE,
                               ndim = 30,
                               gene.list = NULL,
                               hvg.ngenes = 2000,
                               normalization.method = "Log1p",
                               ncores = 1,
                               verbose = TRUE){
   
   slot <- "single-cell"
   # adjust to get processing without removing any cell
   scTypeEval@metadata <- 
      scTypeEval@metadata |>
      dplyr::mutate(k = sample(c("a", "b"),
                               nrow(scTypeEval@metadata),
                               replace = T),
                    sample = .data[[sample]]
      ) |>
      dplyr::select(k, sample)
   
   scTypeEval <- Run.ProcessingData(scTypeEval,
                                    aggregation = slot,
                                    normalization.method = normalization.method,
                                    ident = "k",
                                    sample = "sample",
                                    verbose = verbose,
                                    min.samples = 1,
                                    min.cells = 1)
   
   if(is.null(gene.list)){
      if(verbose){message("No gene list provided, computing HVG")}
      scTypeEval <- Run.HVG(scTypeEval,
                            ngenes = hvg.ngenes,
                            aggregation = slot,
                            ncores = ncores, 
                            verbose = verbose)
   } else {
      if(verbose){message("Using provided gene list: ", names(gene.list))}
      scTypeEval <- add.GeneList(scTypeEval, gene.list = gene.list)
   }
   
   if(reduction){
      if(verbose){message("Producing PCA embeddings")}
      scTypeEval <- Run.PCA(scTypeEval, 
                            ndim = ndim,
                            verbose = F)
      X <- scTypeEval@reductions[[slot]]@embeddings
   } else {
      if(verbose){
         "Not running clustering on low dimensional space, subetting gene list and black list genes if provided."}
      
      mat_ident <- scTypeEval@data[[slot]]
      if(is.null(mat_ident)){
         stop("No processed data slot found for ", slot ,
              ". Please run before `Run.ProcessingData()` or add a data assay.\n")
      }
      gene.list <- .check_genelist(scTypeEval, gene.list, verbose = verbose)
      black.list <- .check_blacklist(scTypeEval, black.list, verbose = verbose)
      mat_ident <- .general_filtering(mat_ident,
                                      black.list = black.list,
                                      gene.list = gene.list,
                                      verbose = verbose)
      X <- mat_ident@matrix
   }
   
   ret <- list(X = X,
               gene.list = scTypeEval@gene.lists)
   return(ret)
}


# auxiliary function to run scTypeEval from processing, dissimilarity to finally obtain consistency

compute_consistency <- function(scTypeEval,
                                ident,
                                sample = "sample",
                                gene.list = NULL,
                                consistency_method,
                                min.samples = 5,
                                min.cells = 10,
                                ncores = 1,
                                verbose = F) {
   
   dissimilarity.method <- sapply(consistency_method, \(x) strsplit(x, " [|] ")[[1]][2])
   consistency.metric  <- sapply(consistency_method, \(x) strsplit(x, " [|] ")[[1]][1])
   
   # all expected clusters
   all_clusters <- purge_label(unique(scTypeEval@metadata[[ident]]))
   all_clusters <- all_clusters[!is.na(all_clusters)]
   
   diss <- wrapper_scTypeEval(
      scTypeEval = scTypeEval,
      aggregation = NULL,
      ident = ident,
      sample = sample,
      gene.list = gene.list,
      dissimilarity.method = dissimilarity.method,
      min.samples = min.samples,
      min.cells = min.cells,
      ncores = ncores,
      verbose = verbose
   )
   
   cons <- get.Consistency(
      diss,
      dissimilarity.slot = dissimilarity.method,
      Consistency.metric = consistency.metric,
      verbose = verbose
   ) |>
      dplyr::mutate(
         consistency.method = paste(consistency.metric,
                                    dissimilarity_method,
                                    sep = " | ")
      ) |> 
      dplyr::select(-consistency.metric, -dissimilarity_method) |>
      dplyr::filter(consistency.method %in% consistency_method) |>
      tidyr::pivot_wider(
         names_from = "consistency.method",
         values_from = "measure"
      )
   
   # product
   if (length(consistency_method) > 1) {
      cons$product <- Reduce(`*`, cons[consistency_method])
   } else {
      cons$product <- cons[[consistency_method]]
   }
   cons
}


# function to split clusters

get.clusters <- function(
      X,
      clustering_method = c("kmeans", "louvain", "leiden"),
      nclusters = 2,
      nstart = 30,
      k = 20,
      transposed = TRUE,
      resolution = 1.0,
      ncores = 1,
      bparam = NULL,
      progressbar = FALSE,
      ...
) {
   clustering_method <- match.arg(clustering_method)
   
   if(clustering_method %in% c("louvain", "leiden")){
      # set paralelization
      param <- set_parallel_params(ncores = ncores,
                                   bparam = bparam,
                                   progressbar = progressbar)
      g <- scran::buildSNNGraph(X,
                                k = k,
                                transposed = transposed,
                                BPPARAM = param)
   }
   
   switch(
      clustering_method,
      
      "kmeans" = {
         return(kmeans(X, centers = nclusters, nstart = nstart)$cluster)
      },
      
      "louvain" = {
         return(cluster_louvain_wrapper(g, resolution = resolution, ...))
      },
      
      "leiden" = {
         return(cluster_leiden_wrapper(g, resolution = resolution, ...))
      },
      
      stop("Unsupported clustering method: ", clustering_method)
   )
}

cluster_louvain_wrapper <- function(g, ...) {
   
   cl <- igraph::cluster_louvain(g, ...)
   labs <- igraph::membership(cl)
   
   return(labs)
}

cluster_leiden_wrapper <- function(g, ...) {
   
   cl <- igraph::cluster_leiden(g, n_iterations = 10, ...)
   labs <- igraph::membership(cl)
   
   return(labs)
}

# function to programmatically set up clustering resolution

pick_resolution <- function(ncells, target = 2) {
   target / sqrt(ncells)
}

