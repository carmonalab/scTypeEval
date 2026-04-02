# auxiliary function to process all dataets, without filtering for obtaining the PCA or 
# matrix with relevant genes on the gene list

process_clustering <- function(scTypeEval,
                               sample = "sample",
                               reduction = TRUE,
                               ndim = 30,
                               gene_list = NULL,
                               hvg_ngenes = 2000,
                               normalization_method = "Log1p",
                               ncores = 1,
                               verbose = TRUE){
   
   slot <- "single-cell"
   # adjust to get processing without removing any cell
   scTypeEval@metadata <- 
      scTypeEval@metadata |>
      dplyr::mutate(k = sample(c("a", "b"),
                               nrow(scTypeEval@metadata),
                               replace = TRUE),
                    sample = .data[[sample]]
      ) |>
      dplyr::select(k, sample)
   
   scTypeEval <- run_processing_data(scTypeEval,
                                    aggregation = slot,
                                    normalization_method = normalization_method,
                                    ident = "k",
                                    sample = "sample",
                                    verbose = verbose,
                                    min_samples = 1,
                                    min_cells = 1)
   
   if(is.null(gene_list)){
      if(verbose){message("No gene list provided, computing HVG")}
      scTypeEval <- run_hvg(scTypeEval,
                            ngenes = hvg_ngenes,
                            aggregation = slot,
                            ncores = ncores, 
                            verbose = verbose)
   } else {
      if(verbose){message("Using provided gene list: ", names(gene_list))}
      scTypeEval <- add_gene_list(scTypeEval, gene_list = gene_list)
   }
   
   if(reduction){
      if(verbose){message("Producing PCA embeddings")}
      scTypeEval <- run_pca(scTypeEval, 
                            ndim = ndim,
                            verbose = FALSE)
      x <- scTypeEval@reductions[[slot]]@embeddings
   } else {
      if(verbose){
         "Not running clustering on low dimensional space, subetting gene list and black list genes if provided."}
      
      mat_ident <- scTypeEval@data[[slot]]
      if(is.null(mat_ident)){
         stop("No processed data slot found for ", slot ,
              ". Please run before `run_processing_data()` or add a data assay.\n")
      }
      gene_list <- check_genelist(scTypeEval, gene_list, verbose = verbose)
      black_list <- check_blacklist(scTypeEval, NULL, verbose = verbose)
      mat_ident <- general_filtering(mat_ident,
                       black_list = black_list,
                       gene_list = gene_list,
                                      verbose = verbose)
      x <- mat_ident@matrix
   }
   
   ret <- list(x = x,
            gene_list = scTypeEval@gene_lists)
   return(ret)
}


# auxiliary function to run scTypeEval from processing, dissimilarity to finally obtain consistency

compute_consistency <- function(scTypeEval,
                                ident,
                                sample = "sample",
                                gene_list = NULL,
                                consistency_method,
                                min_samples = 5,
                                min_cells = 10,
                                ncores = 1,
                                verbose = FALSE) {
   
   dissimilarity_method <- vapply(consistency_method, \(x) strsplit(x, " [|] ")[[1]][2], FUN.VALUE = character(1))
   consistency_metric  <- vapply(consistency_method, \(x) strsplit(x, " [|] ")[[1]][1], FUN.VALUE = character(1))
   
   # all expected clusters
   all_clusters <- purge_label(unique(scTypeEval@metadata[[ident]]))
   all_clusters <- all_clusters[!is.na(all_clusters)]
   
   diss <- wrapper_scTypeEval(
      scTypeEval = scTypeEval,
      aggregation = NULL,
      ident = ident,
      sample = sample,
      gene_list = gene_list,
      dissimilarity_method = dissimilarity_method,
      min_samples = min_samples,
      min_cells = min_cells,
      ncores = ncores,
      verbose = verbose
   )
   
   cons <- get_consistency(
      diss,
      dissimilarity_slot = dissimilarity_method,
      consistency_metric = consistency_metric,
      verbose = verbose
   ) |>
      dplyr::mutate(
         consistency.method = paste(consistency_metric,
                                    dissimilarity_method,
                                    sep = " | ")
      ) |> 
      dplyr::select(-consistency_metric, -dissimilarity_method) |>
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

get_clusters <- function(
   x,
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
      if (!requireNamespace("igraph", quietly = TRUE)) {
         stop("The 'igraph' package is required for ", clustering_method, " clustering. Please install it.",
              call. = FALSE)
      }
      # set paralelization
      param <- set_parallel_params(ncores = ncores,
                                   bparam = bparam,
                                   progressbar = progressbar)
      g <- scran::buildSNNGraph(x,
                                k = k,
                                transposed = transposed,
                                BPPARAM = param)
   }
   
   switch(
      clustering_method,
      
      "kmeans" = {
         return(kmeans(x, centers = nclusters, nstart = nstart)$cluster)
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

