
downsample_factor_level <- function(df, factor_col, level, threshold, seed = 22) {
   set.seed(seed) # Set seed for reproducibility
   
   # Split the dataframe into the target level and the rest
   target_df <- df |> dplyr::filter(!!rlang::sym(factor_col) == level)
   rest_df <- df |> dplyr::filter(!!rlang::sym(factor_col) != level)
   
   # Downsample the target level if it exceeds the threshold
   if (nrow(target_df) > threshold) {
      target_df <- target_df |> 
         dplyr::slice_sample(n = threshold)
   }
   
   # Combine the downsampled target level with the rest of the data
   result_df <- dplyr::bind_rows(target_df, rest_df)
   
   return(result_df)
}



rand.shuffling_group <- function(vector,
                                 group,
                                 rate = 0.1, # proportion of labels to shuffle within each group
                                 seed = 22) {
   set.seed(seed)
   
   nvector <- vector
   
   # Iterate over each unique group
   for (g in unique(group)) {
      idx <- which(group == g)  # indices for the current group
      n <- length(idx)
      
      # Number of labels to shuffle in this group
      shuffle_n <- floor(n * rate)
      
      if (shuffle_n > 0) {
         its <- sample(idx, shuffle_n)  # positions to replace within this group
         nvector[its] <- sample(vector[idx], shuffle_n)  # shuffle labels within group
      }
   }
   
   return(nvector)
}



trim_mean_se <- function(x,
                         trim = 0) {
   n <- length(x)
   # Compute trimmed arithmetic mean
   mean_val <- mean(x, trim = trim)
   
   # Compute standard error
   se <- sd(x) / sqrt(n)
   
   # Return a data frame with mean and SE range
   return(data.frame(
      y = mean_val,
      ymin = mean_val - se,
      ymax = mean_val + se
   ))
}

geo_mean <- function(x) {
   # Replace zeros with a small positive value
   x[x == 0] <- 1e-10
   exp(mean(log(x)))
}

geo_mean_se <- function(x) {
   # Replace zeros with a small positive value
   x[x == 0] <- 1e-10
   n <- length(x)
   
   # Compute geometric mean
   geo_mean <- exp(mean(log(x)))
   
   # Compute standard error
   se <- geo_mean * sd(log(x)) / sqrt(n)
   
   # Return a data frame with mean and SE range
   return(data.frame(
      y = geo_mean,
      ymin = geo_mean - se,
      ymax = geo_mean + se
   ))
}


wrapper_dissimilarity <- function(scTypeEval,
                                  ident,
                                  sample,
                                  gene.list,
                                  reduction = TRUE,
                                  ndim = 30,
                                  normalization.method = "Log1p",
                                  dissimilarity.method = c("WasserStein", "Pseudobulk:Euclidean",
                                                           "Pseudobulk:Cosine", "Pseudobulk:Pearson",
                                                           "BestHit:Match", "BestHit:Score"),
                                  min.samples = 5,
                                  min.cells = 10,
                                  bparam = BiocParallel::SerialParam(),
                                  verbose = TRUE){
   
   scTypeEval <- Run.ProcessingData(scTypeEval,
                                    ident = ident,
                                    sample = sample,
                                    normalization.method = normalization.method,
                                    min.samples = min.samples,
                                    min.cells = min.cells,
                                    verbose = verbose)
   scTypeEval <- add.GeneList(scTypeEval, gene.list)
   
   if(reduction){
      scTypeEval <- Run.PCA(scTypeEval,
                            ndim = ndim)
   }
   
   for(m in dissimilarity.method){
      if(verbose){message(">.  Running ", m, "\n")}
      scTypeEval <- Run.Dissimilarity(scTypeEval,
                                      reduction = reduction,
                                      method = m,
                                      bparam = bparam,
                                      verbose = verbose)
   }
   
   return(scTypeEval)
   
}

wrapper_plots <- function(scTypeEval,
                          reduction.slot = "all", 
                          dissimilarity.slot = "all",
                          reduction = TRUE,
                          label = TRUE,
                          dims = c(1,2),
                          show.legend = FALSE,
                          dir.path,
                          height =15,
                          width = 18,
                          ...){
   if(reduction){
      pca <- plot.PCA(scTypeEval,
                      reduction.slot = reduction.slot,
                      label = label,
                      dims = dims,
                      show.legend = show.legend)
   } else {
      pca <- c()
   }
   
   mds <- plot.MDS(scTypeEval,
                   dissimilarity.slot = dissimilarity.slot,
                   label = label,
                   dims = dims,
                   show.legend = show.legend)
   
   ph <- plot.Heatmap(scTypeEval,
                      dissimilarity.slot = dissimilarity.slot,
                      ...)
   
   # Export as PDF
   dir <- paste0(dir.path, ".pdf")
   # List of plots to print
   plot_list <- c(pca, mds, ph)
   pdf(dir, width = width, height = height)
   
   for (pl in plot_list) {
      if (inherits(pl, "pheatmap")) {
         grid::grid.newpage()
         grid::grid.draw(pl$gtable)
      } else if (inherits(pl, "ggplot")) {
         print(pl)
      } else {
         warning("Unknown plot type: skipping.")
      }
   }
   
   dev.off()
}


# wrapper to evaluate the missclassifiction rate 
wr.missclasify <- function(count_matrix,
                           metadata,
                           ident,
                           sample,
                           rates = c(1, 0.9, 0.75, 0.5, 0.25, 0), # proportion of cell labels kept
                           replicates = 3,
                           gene.list = NULL,
                           reduction = TRUE,
                           ndim = 30,
                           black.list = NULL,
                           dir = NULL,
                           normalization.method = "Log1p",
                           dissimilarity.method = c("WasserStein", "Pseudobulk:Euclidean",
                                                    "Pseudobulk:Cosine", "Pseudobulk:Pearson",
                                                    "BestHit:Match", "BestHit:Score"),
                           IntVal.metric = c("silhouette", "NeighborhoodPurity",
                                             "ward.PropMatch", "Orbital.medoid",
                                             "Average.similarity"),
                           min.samples = 5,
                           min.cells = 10,
                           ncores = 1,
                           bparam = NULL,
                           progressbar = FALSE,
                           seed = 22,
                           KNNGraph_k = 5,
                           hclust.method = "ward.D2",
                           save.plots = TRUE,
                           verbose = TRUE){
   
   if(!is.null(dir)){
      if(verbose){message("\nResults will be stored at ", dir)}
      dir.create(dir, showWarnings = F)
      
      if(save.plots){
         pca.dir <- file.path(dir, "Plots_Missclassify")
         dir.create(pca.dir, showWarnings = F)
      }
   }
   
   # produce missclassification on metadata
   # get seed for each replicate
   sds <- seed + 1:replicates
   names(sds) <- 1:replicates
   
   # vector for annotations
   annotations <- c()
   
   for(s in names(sds)){
      original_vector <- metadata[[ident]]
      
      for(r in rates){
         ra <- 1-r # proportion of cells to shuffle
         new_vector <- rand.shuffling_group(vector = original_vector,
                                            group = metadata[[sample]],
                                            rate = ra,
                                            seed = sds[[s]])
         new_name <- paste("R", s, r, ident, sep = "_")
         metadata[[new_name]] <- new_vector
         annotations <- c(annotations, new_name)
         
      }
   }
   
   ## Create sc object
   sc <- create.scTypeEval(matrix = count_matrix,
                           metadata = metadata,
                           active.ident = ident,
                           black.list = black.list)
   # get the gene list, the same for every run
   if(is.null(gene.list)){
      sc.gl <- Run.ProcessingData(sc,
                                  sample = sample,
                                  normalization.method = normalization.method,
                                  min.samples = min.samples,
                                  min.cells = min.cells,
                                  verbose = verbose)
      sc.gl <- Run.HVG(sc.gl,
                       ncores = ncores)
      gl <- sc.gl@gene.lists
   } else {
      gl <- gene.list
   }
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   if(verbose){message("Running loop of annotations ")}
   df.res <- BiocParallel::bplapply(annotations,
                                    BPPARAM = BiocParallel::SerialParam(),
                                    function(ann){
                                       tryCatch(
                                          {
                                             sc.tmp <- wrapper_dissimilarity(sc,
                                                                             ident = ann,
                                                                             sample = sample,
                                                                             gene.list = gl,
                                                                             reduction = reduction,
                                                                             ndim = ndim,
                                                                             normalization.method = normalization.method,
                                                                             dissimilarity.method = dissimilarity.method,
                                                                             min.samples = min.samples,
                                                                             min.cells = min.cells,
                                                                             bparam = param,
                                                                             verbose = verbose
                                                                             )
                                             # data.frame with consistency outcome
                                             res <- get.Consistency(sc.tmp)
                                             
                                             # accommodate extra data
                                             res <- res |>
                                                dplyr::mutate(rate = as.numeric(as.character(strsplit(ann, "_")[[1]][3])),
                                                              rep = strsplit(ann, "_")[[1]][2],
                                                              original.ident = !!ident,
                                                              task = "Missclassification"
                                                )
                                             
                                             # render PCAs
                                             # only produce for one replicate of the seeds
                                             if(stringr::str_split(ann, "_")[[1]][2] == 1){
                                                # save pdf if indicated
                                                if(save.plots){
                                                   if(verbose){message("\nProducing Plots for ", ann, "\n")}
                                                   fp <- file.path(pca.dir, ann)
                                                   wrapper_plots(sc.tmp,
                                                                 dir.path = fp,
                                                                 reduction = reduction)
                                                }
                                             }
                                             
                                             return(res)
                                          },
                                          error = function(e){
                                             message("\n-X- scTypeEval failed for", ann, "\n")
                                             message(e)
                                          }
                                       )
                                       
                                    })
   
   # concatenate all results
   df.res <- do.call(rbind, df.res)
   
   if(!is.null(dir)){
      saveRDS(df.res,
              file.path(dir,
                        paste0("Missclassification_", ident, ".rds")))
   }
   return(df.res)
   
}

# wrapper to evaluate the number of samples
wr.NSamples <- function(count_matrix,
                        metadata,
                        ident,
                        sample,
                        rates = c(1, 0.9, 0.7, 0.5), # proportion of total samples kept
                        replicates = 3,
                        gene.list = NULL,
                        pca = FALSE,
                        ndim = 30,
                        black.list = NULL,
                        dir = NULL,
                        normalization.method = "Log1p",
                        distance.method = "euclidean",
                        data.type = c("pseudobulk", "pseudobulk_1vsall", "sc"),
                        IntVal.metric = c("silhouette", "NeighborhoodPurity",
                                          "ward.PropMatch", "Leiden.PropMatch",
                                          "modularity"),
                        BH.method = c("Mutual.Score", "Mutual.Match"),
                        classifier = c("SingleR", "Spearman_correlation"),
                        min.samples = 5,
                        min.cells = 10,
                        ncores = 1,
                        bparam = NULL,
                        progressbar = FALSE,
                        seed = 22,
                        k.sc = 20,
                        k.psblk = 5,
                        save.plots = TRUE,
                        dims = c(1,2),
                        show.legend = FALSE,
                        label = TRUE,
                        verbose = TRUE){
   
   if(!is.null(dir)){
      if(verbose){message("\nResults will be stored at ", dir)}
      dir.create(dir, showWarnings = F)
      
      if(save.plots){
         pca.dir <- file.path(dir, "PCA_NSamples")
         dir.create(pca.dir, showWarnings = F)
      }
   }
   
   
   # produce diferents combinations of samples
   # get seed for each replicate
   sds <- seed + 1:replicates
   names(sds) <- 1:replicates
   
   # number of samples
   ss <- unique(metadata[[sample]])
   nsamples <- floor(length(ss) * rates)
   nss <- list()
   
   for(s in names(sds)){
      for(r in nsamples){
         set.seed(sds[[s]])
         dos <- sample(ss, size = r)
         nn <- paste("N", s, r, ident, sep = "_")
         nss[[nn]] <- dos
      }
   }
   
   # gene list should be the same in all the iterations
   if(is.null(gene.list)){
      ## Create sc object
      sc <- create.scTypeEval(matrix = count_matrix,
                              metadata = metadata,
                              black.list = bl)
      sc <- add.HVG(sc,
                    sample = sample,
                    ncores = ncores,
                    progressbar = progressbar)
      gl <- sc@gene.lists
   } else {
      gl <- gene.list
   }
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   df.res <- BiocParallel::bplapply(names(nss),
                                    BPPARAM = param,
                                    function(ns){
                                       tryCatch(
                                          {
                                             ## run scTypeEval
                                             # create scTypeEval object with the number of samples
                                             md <- metadata[metadata[[sample]] %in% nss[[ns]],]
                                             sct <- create.scTypeEval(matrix = count_matrix[,rownames(md)],
                                                                      metadata = md,
                                                                      gene.lists = gl,
                                                                      black.list = bl)
                                             
                                             res <- Run.scTypeEval(scTypeEval = sct,
                                                                   ident = ident,
                                                                   ident_GroundTruth = ident,
                                                                   sample = sample,
                                                                   normalization.method = normalization.method,
                                                                   gene.list = NULL, # fixed gene list earlier
                                                                   pca = pca,
                                                                   ndim = ndim,
                                                                   distance.method = distance.method,
                                                                   IntVal.metric = IntVal.metric,
                                                                   BH.method = BH.method,
                                                                   classifier = classifier,
                                                                   data.type = data.type,
                                                                   min.samples = min.samples,
                                                                   min.cells = min.cells,
                                                                   k.sc = k.sc,
                                                                   k.psblk = k.psblk,
                                                                   black.list = black.list,
                                                                   ncores = 1,
                                                                   bparam = NULL,
                                                                   progressbar = progressbar,
                                                                   verbose = verbose)
                                             
                                             # accommodate extra data
                                             res <- res |>
                                                dplyr::mutate(rate = as.numeric(as.character(strsplit(ns, "_")[[1]][3])),
                                                              rep = strsplit(ns, "_")[[1]][2],
                                                              original.ident = ident,
                                                              task = "NSamples"
                                                )
                                             
                                             # render PCAs
                                             # only produce for one replicate of the seeds
                                             if(strsplit(ns, "_")[[1]][2] == 1){
                                                # save pdf if indicated
                                                if(save.plots){
                                                   if(verbose){message("\nProducing PCAs for ", ns, "\n")}
                                                   pcs <- wr.pca(scTypeEval = sct,
                                                                 ident = ident,
                                                                 sample = sample,
                                                                 normalization.method = normalization.method,
                                                                 gene.list = NULL, # fixed gene list earlier
                                                                 data.type = data.type,
                                                                 min.samples = min.samples,
                                                                 min.cells = min.cells,
                                                                 black.list = black.list,
                                                                 ndim = 50,
                                                                 dims = dims,
                                                                 show.legend = show.legend,
                                                                 label = label,
                                                                 ncores = 1,
                                                                 bparam = NULL,
                                                                 progressbar = progressbar,
                                                                 verbose = verbose)
                                                   
                                                   
                                                   for (p in names(pcs)){
                                                      pdf(file.path(pca.dir, paste0(p, "_", ns, ".pdf")),
                                                          width = 4, height = 4)
                                                      print(pcs[[p]])
                                                      dev.off()
                                                   }
                                                }
                                             }
                                             return(res)
                                          },
                                          error = function(e){
                                             message("\n-X- scTypeEval failed for", ns, "\n")
                                             message(e)
                                          }
                                       )
                                    })
   
   # concatenate all results
   df.res <- do.call(rbind, df.res)
   
   if(!is.null(dir)){
      saveRDS(df.res,
              file.path(dir,
                        paste0("NSamples_", ident, ".rds")))
   }
   
   return(df.res)
   
}

# wrapper to evaluate the consistency metrics when excluding some cell types
wr.Nct <- function(count_matrix,
                   metadata,
                   ident,
                   sample,
                   gene.list = NULL,
                   pca = FALSE,
                   ndim = 30,
                   black.list = NULL,
                   dir = NULL,
                   normalization.method = "Log1p",
                   distance.method = "euclidean",
                   data.type = c("pseudobulk", "pseudobulk_1vsall", "sc"),
                   IntVal.metric = c("silhouette", "NeighborhoodPurity",
                                     "ward.PropMatch", "Leiden.PropMatch",
                                     "modularity"),
                   BH.method = c("Mutual.Score", "Mutual.Match"),
                   classifier = c("SingleR", "Spearman_correlation"),
                   min.samples = 5,
                   min.cells = 10,
                   ncores = 1,
                   bparam = NULL,
                   progressbar = FALSE,
                   k.sc = 20,
                   k.psblk = 5,
                   down.sample.comb = 30,
                   seed = 22,
                   save.plots = TRUE,
                   dims = c(1,2),
                   show.legend = FALSE,
                   label = TRUE,
                   verbose = TRUE){
   
   if(!is.null(dir)){
      if(verbose){message("Results will be stored at ", dir)}
      dir.create(dir, showWarnings = F)
      
      if(save.plots){
         pca.dir <- file.path(dir, "PCA_Nct")
         dir.create(pca.dir, showWarnings = F)
      }
   }
   
   
   # get all possible combinations
   allcts <- unique(metadata[[ident]])
   allcts <- allcts[!is.na(allcts)]
   
   
   # Generate combinations of lengths 2 to n-1
   cts <- sample_variable_length_combinations(allcts,
                                              num_samples = down.sample.comb,
                                              seed = seed)
   # add all cell types
   allelements <- paste(allcts, collapse = "-")
   cts[[allelements]] <- allcts
   
   # gene list should be the same in all the iterations
   if(is.null(gene.list)){
      ## Create sc object
      sc <- create.scTypeEval(matrix = count_matrix,
                              metadata = metadata,
                              black.list = bl)
      sc <- add.HVG(sc,
                    sample = sample,
                    ncores = ncores,
                    progressbar = progressbar)
      gl <- sc@gene.lists
   } else {
      gl <- gene.list
   }
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   df.res <- BiocParallel::bplapply(names(cts),
                                    BPPARAM = param,
                                    function(ns){
                                       tryCatch(
                                          {
                                             # create scTypeEval object with the number of samples
                                             md <- metadata[metadata[[ident]] %in% cts[[ns]],]
                                             sct <- create.scTypeEval(matrix = count_matrix[,rownames(md)],
                                                                      metadata = md,
                                                                      gene.lists = gl,
                                                                      black.list = bl)
                                             
                                             res <- Run.scTypeEval(scTypeEval = sct,
                                                                   ident = ident,
                                                                   ident_GroundTruth = ident,
                                                                   sample = sample,
                                                                   normalization.method = normalization.method,
                                                                   gene.list = NULL, # fixed gene list earlier
                                                                   pca = pca,
                                                                   ndim = ndim,
                                                                   distance.method = distance.method,
                                                                   IntVal.metric = IntVal.metric,
                                                                   BH.method = BH.method,
                                                                   classifier = classifier,
                                                                   data.type = data.type,
                                                                   min.samples = min.samples,
                                                                   min.cells = min.cells,
                                                                   k.sc = k.sc,
                                                                   k.psblk = k.psblk,
                                                                   black.list = black.list,
                                                                   ncores = 1,
                                                                   bparam = NULL,
                                                                   progressbar = progressbar,
                                                                   verbose = verbose)
                                             
                                             # accommodate extra data
                                             res <- res |>
                                                dplyr::mutate(rate = ns,
                                                              rep = NA,
                                                              original.ident = ident,
                                                              task = "NSamples"
                                                )
                                             
                                             # render PCAs
                                             if(save.plots){
                                                if(verbose){message("\nProducing PCAs for ", ns, "\n")}
                                                pcs <- wr.pca(scTypeEval = sct,
                                                              ident = ident,
                                                              sample = sample,
                                                              normalization.method = normalization.method,
                                                              gene.list = NULL, # fixed gene list earlier
                                                              data.type = data.type,
                                                              min.samples = min.samples,
                                                              min.cells = min.cells,
                                                              black.list = black.list,
                                                              ndim = 50,
                                                              dims = dims,
                                                              show.legend = show.legend,
                                                              label = label,
                                                              ncores = 1,
                                                              bparam = NULL,
                                                              progressbar = progressbar,
                                                              verbose = verbose)
                                                
                                                for (p in names(pcs)){
                                                   pdf(file.path(pca.dir, paste0(p, "_", ns, ".pdf")),
                                                       width = 4, height = 4)
                                                   print(pcs[[p]])
                                                   dev.off()
                                                }
                                                
                                             }
                                             return(res)
                                          },
                                          error = function(e){
                                             message("\n-X- scTypeEval failed for", ns, "\n")
                                             message(e)
                                          }
                                       )
                                       
                                    })
   
   # concatenate all results
   df.res <- do.call(rbind, df.res) |>
      # convert to factor to then compute the fit.constant
      dplyr::mutate(rate = factor(rate,
                                  levels = names(cts))
      )
   
   if(!is.null(dir)){
      saveRDS(df.res,
              file.path(dir,
                        paste0("Nct_", ident, ".rds")))
   }
   
   return(df.res)
}

# wrapper to evaluate the number of cells in a cell type
wr.NCell <- function(count_matrix,
                     metadata,
                     ident,
                     sample,
                     rates = c(1, 0.75, 0.5, 0.25), # proportion of total samples kept
                     ctype = NULL, # which cell type lower the proportion
                     replicates = 3,
                     gene.list = NULL,
                     pca = FALSE,
                     ndim = 30,
                     black.list = NULL,
                     dir = NULL,
                     normalization.method = "Log1p",
                     distance.method = "euclidean",
                     data.type = c("pseudobulk", "pseudobulk_1vsall", "sc"),
                     IntVal.metric = c("silhouette", "NeighborhoodPurity",
                                       "ward.PropMatch", "Leiden.PropMatch",
                                       "modularity"),
                     BH.method = c("Mutual.Score", "Mutual.Match"),
                     classifier = c("SingleR", "Spearman_correlation"),
                     min.samples = 5,
                     min.cells = 10,
                     ncores = 1,
                     bparam = NULL,
                     progressbar = FALSE,
                     seed = 22,
                     k.sc = 20,
                     k.psblk = 5,
                     save.plots = TRUE,
                     dims = c(1,2),
                     show.legend = FALSE,
                     label = TRUE,
                     verbose = TRUE){
   
   if(!is.null(dir)){
      if(verbose){message("\nResults will be stored at ", dir)}
      dir.create(dir, showWarnings = F)
      
      if(save.plots){
         pca.dir <- file.path(dir, "PCA_NCell")
         dir.create(pca.dir, showWarnings = F)
      }
   }
   
   if(is.null(ctype)){
      ctype <- unique(metadata[[ident]])[1]
   }
   if(verbose){message("\nReducing number of cells for ", ctype, "\n")}
   
   
   # produce diferents combinations of samples
   # get seed for each replicate
   sds <- seed + 1:replicates
   names(sds) <- 1:replicates
   
   # get number of cells from ctype
   nb <- table(metadata[[ident]])
   nb <- nb[ctype]
   
   
   # gene list should be the same in all the iterations
   if(is.null(gene.list)){
      ## Create sc object
      sc <- create.scTypeEval(matrix = count_matrix,
                              metadata = metadata,
                              black.list = bl)
      sc <- add.HVG(sc,
                    sample = sample,
                    ncores = ncores,
                    progressbar = progressbar)
      gl <- sc@gene.lists
   } else {
      gl <- gene.list
   }
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   combi <- expand.grid(names(sds), rates)
   
   df.res <- BiocParallel::bplapply(1:nrow(combi),
                                    BPPARAM = param,
                                    function(i){
                                       tryCatch(
                                          {
                                             ## run scTypeEval
                                             s <- combi[i,][1] |> as.numeric()
                                             ns <- combi[i,][2] |> as.numeric()
                                             # create scTypeEval object with the number of samples
                                             # subset cell type
                                             th <- floor(ns*nb)
                                             # subset metadata for specific cell
                                             md <- downsample_factor_level(metadata,
                                                                           ident,
                                                                           ctype,
                                                                           threshold = th,
                                                                           seed = sds[[s]])
                                             #create scTypeEval object
                                             sct <- create.scTypeEval(matrix = count_matrix[,rownames(md)],
                                                                      metadata = md,
                                                                      gene.lists = gl,
                                                                      black.list = bl)
                                             
                                             res <- Run.scTypeEval(scTypeEval = sct,
                                                                   ident = ident,
                                                                   ident_GroundTruth = ident,
                                                                   sample = sample,
                                                                   normalization.method = normalization.method,
                                                                   gene.list = NULL, # fixed gene list earlier
                                                                   pca = pca,
                                                                   ndim = ndim,
                                                                   distance.method = distance.method,
                                                                   IntVal.metric = IntVal.metric,
                                                                   BH.method = BH.method,
                                                                   classifier = classifier,
                                                                   data.type = data.type,
                                                                   min.samples = min.samples,
                                                                   min.cells = min.cells,
                                                                   k.sc = k.sc,
                                                                   k.psblk = k.psblk,
                                                                   black.list = black.list,
                                                                   ncores = 1,
                                                                   bparam = NULL,
                                                                   progressbar = progressbar,
                                                                   verbose = verbose)
                                             
                                             # accommodate extra data
                                             res <- res |>
                                                dplyr::mutate(rate = as.numeric(as.character(ns)),
                                                              rep = s,
                                                              original.ident = ident,
                                                              task = "NCell"
                                                )
                                             
                                             # render PCAs
                                             # only produce for one replicate of the seeds
                                             if(s == 1){
                                                # save pdf if indicated
                                                if(save.plots){
                                                   if(verbose){message("\nProducing PCAs for ", as.character(ns), "\n")}
                                                   pcs <- wr.pca(scTypeEval = sct,
                                                                 ident = ident,
                                                                 sample = sample,
                                                                 normalization.method = normalization.method,
                                                                 gene.list = NULL, # fixed gene list earlier
                                                                 data.type = data.type,
                                                                 min.samples = min.samples,
                                                                 min.cells = min.cells,
                                                                 black.list = black.list,
                                                                 ndim = 50,
                                                                 dims = dims,
                                                                 show.legend = show.legend,
                                                                 label = label,
                                                                 ncores = 1,
                                                                 bparam = NULL,
                                                                 progressbar = progressbar,
                                                                 verbose = verbose)
                                                   
                                                   
                                                   for (p in names(pcs)){
                                                      pdf(file.path(pca.dir, paste0(p, "_", as.character(ns), ".pdf")),
                                                          width = 4, height = 4)
                                                      print(pcs[[p]])
                                                      dev.off()
                                                   }
                                                }
                                             }
                                             return(res)
                                          },
                                          error = function(e){
                                             message("\n-X- scTypeEval failed for", i, "\n")
                                             message(e)
                                          }
                                       )
                                    })
   
   # concatenate all results
   df.res <- do.call(rbind, df.res)
   
   if(!is.null(dir)){
      saveRDS(df.res,
              file.path(dir,
                        paste0("NCell_", ctype, "_", ident, ".rds")))
   }
   
   
   return(df.res)
   
}



wr.pca <- function(scTypeEval,
                   ident = NULL,
                   sample = NULL,
                   normalization.method = c("Log1p", "CLR", "pearson"),
                   gene.list = NULL,
                   data.type = c("sc", "pseudobulk",
                                 "pseudobulk_1vsall"),
                   min.samples = 2,
                   min.cells = 10,
                   black.list = NULL,
                   ndim = 50,
                   dims = c(1,2),
                   show.legend = FALSE,
                   label = TRUE,
                   ncores = 1,
                   bparam = NULL,
                   progressbar = FALSE,
                   verbose = TRUE){
   
   pc.list <- c()
   for(dt in data.type){
      tit <- paste(dt, ident, sep = "_")
      
      if(dt == "sc"){spl = NULL } else { spl = sample}
      
      sc.tmp <- add.PCA(scTypeEval,
                        ident = ident,
                        sample = spl,
                        normalization.method = normalization.method,
                        gene.list = gene.list,
                        data.type = dt,
                        min.samples = min.samples,
                        min.cells = min.cells,
                        black.list = black.list,
                        ndim = ndim,
                        ncores = ncores,
                        bparam = bparam,
                        progressbar = progressbar,
                        verbose = TRUE)
      
      pc <- plot.PCA(sc.tmp,
                     dims = dims,
                     show.legend = show.legend,
                     label = label)
      
      pc <- lapply(pc, function(p){
         p + labs(color = ident)
      })
      
      pc.list <- c(pc.list, pc)
      
   }
   
   return(pc.list)
}

# wrapper to hirarchally merge cell types
wr.mergeCT <- function(count_matrix,
                       metadata,
                       ident,
                       sample,
                       gene.list = NULL,
                       pca = FALSE,
                       ndim = 30,
                       black.list = NULL,
                       dir = NULL,
                       normalization.method = "Log1p",
                       distance.method = "euclidean",
                       data.type = c("pseudobulk", "pseudobulk_1vsall", "sc"),
                       IntVal.metric = c("silhouette", "NeighborhoodPurity",
                                         "ward.PropMatch", "Leiden.PropMatch",
                                         "modularity"),
                       BH.method = c("Mutual.Score", "Mutual.Match"),
                       classifier = c("SingleR", "Spearman_correlation"),
                       min.samples = 5,
                       min.cells = 10,
                       ncores = 1,
                       bparam = NULL,
                       progressbar = FALSE,
                       k.sc = 20,
                       k.psblk = 5,
                       save.plots = TRUE,
                       dims = c(1,2),
                       show.legend = FALSE,
                       label = TRUE,
                       verbose = TRUE){
   
   if(!is.null(dir)){
      if(verbose){message("\nResults will be stored at ", dir)}
      dir.create(dir, showWarnings = F)
      
      if(save.plots){
         pca.dir <- file.path(dir, "PCA_mergeCT")
         dir.create(pca.dir, showWarnings = F)
      }
   }
   
   ## Create sc object
   sc <- create.scTypeEval(matrix = count_matrix,
                           metadata = metadata,
                           black.list = black.list)
   
   sc <- add.HVG(sc,
                 sample = sample,
                 ncores = ncores,
                 progressbar = progressbar)
   gl <- sc@gene.lists
   
   original_vector <- factor(gsub(" |_", ".", metadata[[ident]]))
   j <- length(unique(original_vector))
   new_name <- paste("J", j, ident, sep = "_")
   metadata[[new_name]] <- original_vector
   annotations <- new_name
   
   keep <- rownames(sc@counts) %in% gl[["HVG"]]
   mat <- sc@counts[keep,]
   
   # remove black list genes
   mat <- mat[!rownames(mat) %in% sc@black.list,]
   
   mats <- get.matrix(mat,
                      data.type = "sc",
                      ident = original_vector,
                      sample = NULL,
                      min.samples = min.samples,
                      min.cells = min.cells)
   
   # get distance
   norm.matrix <- Normalize_data(mats@matrix,
                                 normalization.method = normalization.method)
   
   idents <- as.character(mats@ident)
   # keep only celltypes fulfulling min.samples and min.cells conditions
   metadata <- metadata[metadata[[new_name]] %in% unique(idents),]
   
   while(length(unique(idents))>2){
      centroids <- compute_centroids(norm.matrix, idents)
      dist <- get.distance(norm.mat = centroids,
                           distance.method = distance.method)
      dist_mat <- as.matrix(dist)
      diag(dist_mat) <- Inf  # Ignore self-distance
      # Find the closest pair of cell types
      min_index <- which(dist_mat == min(dist_mat), arr.ind = TRUE)[1, ]
      type1 <- colnames(centroids)[min_index[1]]
      type2 <- colnames(centroids)[min_index[2]]
      # Merge the two cell types
      merged_name <- paste0(type1, "+", type2)
      idents[idents %in% c(type1, type2)] <- merged_name
      j <- j - 1
      new_name <- paste("J", j, ident, sep = "_")
      metadata[[new_name]] <- idents
      annotations <- c(annotations, new_name)
      
   }
   
   ## ReCreate sc object
   sc <- create.scTypeEval(matrix = count_matrix[,rownames(metadata)],
                           metadata = metadata,
                           black.list = black.list)
   if(is.null(gene.list)){
      sc <- add.HVG(sc,
                    sample = sample,
                    ncores = ncores,
                    progressbar = progressbar)
      gl <- sc@gene.lists
   } else {
      sc@gene.lists <- gene.list
   }
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   df.res <- BiocParallel::bplapply(annotations,
                                    BPPARAM = param,
                                    function(ann){
                                       tryCatch(
                                          {
                                             res <- Run.scTypeEval(scTypeEval = sc,
                                                                   ident = ann,
                                                                   sample = sample,
                                                                   normalization.method = normalization.method,
                                                                   gene.list = NULL, # Run all gene.lists in the object
                                                                   pca = pca,
                                                                   ndim = ndim,
                                                                   distance.method = distance.method,
                                                                   IntVal.metric = IntVal.metric,
                                                                   BH.method = BH.method,
                                                                   classifier = classifier,
                                                                   data.type = data.type,
                                                                   min.samples = min.samples,
                                                                   min.cells = min.cells,
                                                                   k.sc = k.sc,
                                                                   k.psblk = k.psblk,
                                                                   black.list = black.list,
                                                                   ncores = 1,
                                                                   bparam = NULL,
                                                                   progressbar = progressbar,
                                                                   verbose = verbose)
                                             
                                             # accommodate extra data
                                             res <- res |>
                                                dplyr::mutate(
                                                   rate = as.numeric(as.character(strsplit(ann, "_")[[1]][2])),
                                                   rep = NA,
                                                   original.ident = ident,
                                                   task = "mergeCT"
                                                )
                                             
                                             # render PCAs
                                             # save pdf if indicated
                                             if(save.plots){
                                                if(verbose){message("\nProducing PCAs for ", ann, "\n")}
                                                pcs <- wr.pca(scTypeEval = sc,
                                                              ident = ann,
                                                              sample = sample,
                                                              normalization.method = normalization.method,
                                                              gene.list = NULL,
                                                              data.type = data.type,
                                                              min.samples = min.samples,
                                                              min.cells = min.cells,
                                                              black.list = black.list,
                                                              ndim = 50,
                                                              dims = dims,
                                                              show.legend = show.legend,
                                                              label = label,
                                                              ncores = 1,
                                                              bparam = NULL,
                                                              progressbar = progressbar,
                                                              verbose = verbose)
                                                
                                                
                                                for (p in names(pcs)){
                                                   pdf(file.path(pca.dir, paste0(p, "_", ann, ".pdf")),
                                                       width = 4, height = 4)
                                                   print(pcs[[p]])
                                                   dev.off()
                                                }
                                                
                                             }
                                             return(res)
                                          },
                                          error = function(e){
                                             message("\n-X- scTypeEval failed for", ann, "\n")
                                             message(e)
                                          }
                                       )
                                    })
   
   # concatenate all results
   df.res <- do.call(rbind, df.res)
   
   if(!is.null(dir)){
      saveRDS(df.res,
              file.path(dir,
                        paste0("mergeCT_", ident, ".rds")))
   }
   
   return(df.res)
   
}


wr.assayPlot <- function(df,
                         type = c(1,2),
                         return.df = FALSE,
                         verbose = TRUE,
                         x.label = "rate",
                         trim = 0,
                         title = "",
                         combine = TRUE){
   
   rsq <- df |> 
      dplyr::group_by(gene.list, data.type, consistency.metric, rate) |>
      dplyr::summarize(mm = mean(scaled_measure,
                                 trim = trim,
                                 na.rm = T)) |>
      dplyr::group_by(gene.list, data.type, consistency.metric)
   
   # remove NA or NaN
   rsq <- rsq[complete.cases(rsq),]
   
   type <- type[1]
   
   if(type == 1){
      rsq <- rsq |>
         dplyr::summarize(r.squared = round(fit.ReferenceLine(x = rate, y = mm)$r.squared, 3),
                          pval = round(fit.ReferenceLine(x = rate, y = mm)$p.value, 3),
                          p.val_fill = ifelse(pval < 0.05, "sig", "ns")
         )
      type.verb <- "fit.ReferenceLine, expecting a perfect curve of intercept = 0 and slope = 1\n"
      
   } else if (type == 2){
      rsq <- rsq |>
         dplyr::summarize(r.squared = round(fit.Constant(x = as.numeric(rate), y = mm)$`1-rss`, 3),
                          pval = round(fit.Constant(x = as.numeric(rate), y = mm)$p.value, 3),
                          p.val_fill = ifelse(pval < 0.05, "sig", "ns")
         )
      type.verb <- "fit.Constant, expecting a flat curve with slope = 0\n"
   }
   
   if(verbose){message("\nEvaluating metric using ", type.verb)}
   
   if(!return.df){
      dts <- unique(df[["data.type"]])
      pls <- list()
      
      range.x <- c(min(as.numeric(df[["rate"]])),
                   max(as.numeric(df[["rate"]]))
      )
      
      for(d in dts){
         dftmp <- df |>
            dplyr::filter(data.type == d)
         rsqtmp <- rsq |>
            dplyr::filter(data.type == d)
         
         ncol <- length(unique(dftmp[["consistency.metric"]]))
         
         plc <- dftmp |>
            ggplot2::ggplot(ggplot2::aes(rate, scaled_measure)) +
            ggplot2::stat_summary(ggplot2::aes(group = celltype,
                                               color = celltype),
                                  geom = "point",
                                  alpha = 0.6,
                                  shape = 1,
                                  fun = mean) +
            ggplot2::stat_summary(ggplot2::aes(group = consistency.metric),
                                  geom = "line",
                                  fun = function(y){mean(y, trim = trim)}) +
            ggplot2::stat_summary(ggplot2::aes(group = consistency.metric),
                                  geom = "point",
                                  fun = function(y){mean(y, trim = trim)}) +
            ggplot2::geom_label(data = rsqtmp,
                                ggplot2::aes(x = range.x[1] + 0.2, y = 1.2,
                                             label = r.squared,
                                             fill = r.squared),
                                alpha = 0.4,
                                show.legend = FALSE) +
            ggplot2::geom_label(data = rsqtmp,
                                ggplot2::aes(x = (range.x[2] - 0.2), y = 1.2,
                                             label = pval),
                                fill = "white",
                                alpha = 0.4,
                                show.legend = FALSE) + 
            ggplot2::scale_fill_gradient(low = "red", high = "chartreuse3",
                                         breaks = c(-100, -50, 0, 0.25, 0.5, 0.50, 0.75, 1)) +
            ggplot2::geom_hline(yintercept = 1,
                                linetype = "dotted") +
            ggplot2::facet_wrap(~consistency.metric,
                                ncol =  ncol,
                                drop = FALSE) +
            ggplot2::labs(x = x.label,
                          title = d) +
            ggplot2::ylim(c(0,1.3)) +
            ggpubr::theme_classic2() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust =1))
         
         pls[[d]] <- plc
      }
      
      if(combine){
         plg <- ggpubr::ggarrange(plotlist = pls, ncol = 1, common.legend = TRUE)
         plg <- ggpubr::annotate_figure(plg, top = ggpubr::text_grob(title, 
                                                                     color = "darkblue",
                                                                     face = "bold",
                                                                     size = 20))
      } else {
         plg <- pls
      }
      
      return(plg)
   } else {
      return(as.data.frame(rsq))
   }
}

