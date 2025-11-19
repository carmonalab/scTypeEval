
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

rand.Split_group <- function(vector,
                             group,
                             rate = 0.5, # proportion of labels to shuffle within each group
                             seed = 22,
                             celltype = "ctA") {

   
   if(!celltype %in% vector){
      stop("No cell type found in vector")
   }
   
   set.seed(seed)
   
   # identify cells of the chosen type
   idx_ct <- which(vector == celltype)
   
   # working copies
   new_labels <- vector[idx_ct]  # keep original names = "ctA" (spl1)
   sample_sub <- group[idx_ct]   
   
   # iterate over samples
   for (s in unique(sample_sub)) {
      idx_s     <- which(sample_sub == s)
      n_s       <- length(idx_s)
      n_shuffle <- floor(n_s * rate)
      
      if(n_shuffle > 0){
         # choose which cells go to spl2
         chosen <- sample(idx_s, n_shuffle)
         new_labels[chosen] <- paste0(celltype, ".Nspl2")
      }
   }
   
   # assign back
   vector[idx_ct] <- new_labels
   
   return(vector)
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
                                             "Average.similarity", "2label.silhouette"),
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
   
   original_vector <- metadata[[ident]]
   # vector for annotations
   annotations <- c()
   
   for(s in names(sds)){
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
   df.res <- lapply(annotations,
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
                             message("\n-X- scTypeEval failed for ", ann, "\n", e)
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
                        rates = c(1, 0.9, 0.7, 0.5), # proportion of cell labels kept
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
                                          "Average.similarity", "2label.silhouette"),
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
         pca.dir <- file.path(dir, "Plots_Nsamples")
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
         nn <- paste("N", s, r, sep = "_")
         nss[[nn]] <- dos
      }
   }
   
   
   # get the gene list, the same for every run
   if(is.null(gene.list)){
      ## Create sc object
      sc <- create.scTypeEval(matrix = count_matrix,
                              metadata = metadata,
                              active.ident = ident,
                              black.list = black.list)
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
   df.res <- lapply(names(nss),
                    function(ns){
                       tryCatch(
                          {
                             ## run scTypeEval
                             # create scTypeEval object with the number of samples
                             md <- metadata[metadata[[sample]] %in% nss[[ns]],]
                             sc.tmp <- create.scTypeEval(matrix = count_matrix[,rownames(md)],
                                                         metadata = md,
                                                         gene.lists = gl,
                                                         black.list = black.list)
                             
                             sc.tmp <- wrapper_dissimilarity(sc.tmp,
                                                             ident = ident,
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
                                dplyr::mutate(rate = as.numeric(as.character(strsplit(ns, "_")[[1]][3])),
                                              rep = strsplit(ns, "_")[[1]][2],
                                              original.ident = !!ident,
                                              task = "NSamples"
                                )
                             
                             # render PCAs
                             # only produce for one replicate of the seeds
                             if(stringr::str_split(ns, "_")[[1]][2] == 1){
                                # save pdf if indicated
                                if(save.plots){
                                   if(verbose){message("\nProducing Plots for ", ns, "\n")}
                                   fp <- file.path(pca.dir, ns)
                                   wrapper_plots(sc.tmp,
                                                 dir.path = fp,
                                                 reduction = reduction)
                                }
                             }
                             
                             return(res)
                          },
                          error = function(e){
                             message("\n-X- scTypeEval failed for ", ns, "\n", e)
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
                                     "Average.similarity", "2label.silhouette"),
                   min.samples = 5,
                   min.cells = 10,
                   ncores = 1,
                   bparam = NULL,
                   progressbar = FALSE,
                   seed = 22,
                   KNNGraph_k = 5,
                   hclust.method = "ward.D2",
                   save.plots = TRUE,
                   verbose = TRUE,
                   down.sample.comb = 30){
   
   if(!is.null(dir)){
      if(verbose){message("Results will be stored at ", dir)}
      dir.create(dir, showWarnings = F)
      
      if(save.plots){
         pca.dir <- file.path(dir, "Plots_Nct")
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
   
   # get the gene list, the same for every run
   if(is.null(gene.list)){
      ## Create sc object
      sc <- create.scTypeEval(matrix = count_matrix,
                              metadata = metadata,
                              active.ident = ident,
                              black.list = black.list)
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
   
   df.res <- lapply(names(cts),
                    function(ns){
                       tryCatch(
                          {
                             # create scTypeEval object with the number of samples
                             md <- metadata[metadata[[ident]] %in% cts[[ns]],]
                             sc.tmp <- create.scTypeEval(matrix = count_matrix[,rownames(md)],
                                                         metadata = md,
                                                         gene.lists = gl,
                                                         black.list = black.list)
                             sc.tmp <- wrapper_dissimilarity(sc.tmp,
                                                             ident = ident,
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
                                dplyr::mutate(rate = ns,
                                              rep = NA,
                                              original.ident = !!ident,
                                              task = "Nct"
                                )
                             
                             # render plots
                             # save pdf if indicated
                             if(save.plots){
                                if(verbose){message("\nProducing Plots for ", ns, "\n")}
                                fp <- file.path(pca.dir, ns)
                                wrapper_plots(sc.tmp,
                                              dir.path = fp,
                                              reduction = reduction)
                             }
                             
                             return(res)
                          },
                          error = function(e){
                             message("\n-X- scTypeEval failed for ", ns, "\n", e)
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
                                       "Average.similarity", "2label.silhouette"),
                     min.samples = 5,
                     min.cells = 10,
                     ncores = 1,
                     bparam = NULL,
                     progressbar = FALSE,
                     seed = 22,
                     KNNGraph_k = 5,
                     hclust.method = "ward.D2",
                     save.plots = TRUE,
                     verbose = TRUE
){
   
   if(!is.null(dir)){
      if(verbose){message("\nResults will be stored at ", dir)}
      dir.create(dir, showWarnings = F)
      
      if(save.plots){
         pca.dir <- file.path(dir, "Plots_NCell")
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
   
   
   # get the gene list, the same for every run
   if(is.null(gene.list)){
      ## Create sc object
      sc <- create.scTypeEval(matrix = count_matrix,
                              metadata = metadata,
                              active.ident = ident,
                              black.list = black.list)
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
   
   combi <- expand.grid(names(sds), rates)
   
   df.res <- lapply(1:nrow(combi),
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
                             sc.tmp <- create.scTypeEval(matrix = count_matrix[,rownames(md)],
                                                         metadata = md,
                                                         gene.lists = gl,
                                                         black.list = black.list)
                             
                             sc.tmp <- wrapper_dissimilarity(sc.tmp,
                                                             ident = ident,
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
                                dplyr::mutate(rate = as.numeric(as.character(ns)),
                                              rep = s,
                                              original.ident = !!ident,
                                              task = "NCell"
                                )
                             
                             # render PCAs
                             # only produce for one replicate of the seeds
                             if(s == 1){
                                # save pdf if indicated
                                # save pdf if indicated
                                if(save.plots){
                                   if(verbose){message("\nProducing Plots for ", ns, "\n")}
                                   fp <- file.path(pca.dir, ns)
                                   wrapper_plots(sc.tmp,
                                                 dir.path = fp,
                                                 reduction = reduction)
                                }
                             }
                             return(res)
                          },
                          error = function(e){
                             message("\n-X- scTypeEval failed for", i, "\n", e)
                          }
                       )
                    })
   
   # concatenate all results
   df.res <- do.call(rbind, df.res)
   ctype_name <- purge_label(ctype)
   
   if(!is.null(dir)){
      saveRDS(df.res,
              file.path(dir,
                        paste0("NCell_", ctype_name, "_", ident, ".rds")))
   }
   
   
   return(df.res)
   
}



# wrapper to hirarchally merge cell types
wr.mergeCT <- function(count_matrix,
                       metadata,
                       ident,
                       sample,
                       distance.method = "euclidean",
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
                                         "Average.similarity", "2label.silhouette"),
                       min.samples = 5,
                       min.cells = 10,
                       ncores = 1,
                       bparam = NULL,
                       progressbar = FALSE,
                       KNNGraph_k = 5,
                       hclust.method = "ward.D2",
                       save.plots = TRUE,
                       verbose = TRUE){
   
   if(!is.null(dir)){
      if(verbose){message("\nResults will be stored at ", dir)}
      dir.create(dir, showWarnings = F)
      
      if(save.plots){
         pca.dir <- file.path(dir, "Plots_mergeCT")
         dir.create(pca.dir, showWarnings = F)
      }
   }
   
   ## Create sc object
   sc <- create.scTypeEval(matrix = count_matrix,
                           metadata = metadata,
                           active.ident = ident,
                           black.list = black.list)
   
   sc <- Run.ProcessingData(sc,
                            sample = sample,
                            normalization.method = normalization.method,
                            min.samples = min.samples,
                            min.cells = min.cells,
                            verbose = verbose)
   if(is.null(gene.list)){
      
      sc <- Run.HVG(sc,
                    ncores = ncores)
      gl <- sc@gene.lists
   } else {
      gl <- gene.list
   }
   
   gene.list <- .check_genelist(sc, gene.list, verbose = verbose)
   black.list <- .check_blacklist(sc, black.list, verbose = verbose)
   mat_ident <- .general_filtering(sc@data[["single-cell"]],
                                   black.list = black.list,
                                   gene.list = gene.list,
                                   verbose = verbose)
   mat <- mat_ident@matrix
   
   idents <- as.character(mat_ident@ident)
   # keep only celltypes fulfulling min.samples and min.cells conditions
   metadata <- metadata[colnames(mat),]
   
   original_vector <- factor(purge_label(metadata[[ident]]))
   j <- length(unique(original_vector))
   new_name <- paste("J", j, ident, sep = "_")
   metadata[[new_name]] <- original_vector
   annotations <- new_name
   
   while(length(unique(idents))>2){
      centroids <- compute_centroids(norm.mat = mat,
                                     ident = idents)
      dist <- get.distance(norm.mat = centroids,
                           distance.method = distance.method,
                           verbose = F)
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
                           gene.lists = gl,
                           black.list = black.list)
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   df.res <- lapply(annotations,
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
                                dplyr::mutate(
                                   rate = as.numeric(as.character(strsplit(ann, "_")[[1]][2])),
                                   rep = NA,
                                   original.ident = !!ident,
                                   task = "mergeCT"
                                )
                             
                             # render PCAs
                             # save pdf if indicated
                             if(save.plots){
                                if(verbose){message("\nProducing Plots for ", ann, "\n")}
                                fp <- file.path(pca.dir, ann)
                                wrapper_plots(sc.tmp,
                                              dir.path = fp,
                                              reduction = reduction)
                                
                             }
                             return(res)
                          },
                          error = function(e){
                             message("\n-X- scTypeEval failed for", ann, "\n, e")
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

# wrapper to split a cell type into 2 different, but highly similar transcriptomically
wr.splitCellType <- function(count_matrix,
                             metadata,
                             ident,
                             sample,
                             rates = c(1, 0.5),
                             ctype = NULL, # which cell type to split
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
                                               "Average.similarity", "2label.silhouette"),
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
         pca.dir <- file.path(dir, "Plots_SplitCellType")
         dir.create(pca.dir, showWarnings = F)
      }
   }
   
   # produce missclassification on metadata
   # get seed for each replicate
   sds <- seed + 1:replicates
   names(sds) <- 1:replicates
   
   original_vector <- metadata[[ident]]
   
   # get cell type to split, defined or most abundant
   if(is.null(ctype)){
      ctype <- original_vector[which.max(table(original_vector))]
   }
   if(verbose){message("\nReducing number of cells for ", ctype, "\n")}
   
   
   # vector for annotations
   annotations <- c()
   
   for(s in names(sds)){
      for(r in rates){
         ra <- 1-r # proportion of cells to shuffle
         new_vector <- rand.Split_group(vector = original_vector,
                                        group = metadata[[sample]],
                                        rate = ra,
                                        seed = sds[[s]],
                                        celltype = ctype)
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
   df.res <- lapply(annotations,
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
                                              task = "SplitCelltype"
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
                             message("\n-X- scTypeEval failed for ", ann, "\n", e)
                          }
                       )
                       
                    })
   
   # concatenate all results
   df.res <- do.call(rbind, df.res)
   ctype_name <- purge_label(ctype)
   
   if(!is.null(dir)){
      saveRDS(df.res,
              file.path(dir,
                        paste0("SplitCelltype_", ctype_name, "_", ident, ".rds")))
   }
   
   return(df.res)
}


wr.assayPlot <- function(df,
                         type = c("Monotonic", "ReferenceLine", "Constant"),
                         group = c("celltype", "ident"),
                         return.df = FALSE,
                         verbose = TRUE,
                         xlabel = "rate",
                         title = "",
                         combine = TRUE){
   
   group <- group[1] |> tolower()
   
   if(group == "celltype"){
      rsq <- df |> 
         dplyr::group_by(dissimilarity_method, consistency.metric, celltype, rate) |>
         dplyr::summarize(measure = mean(measure)) |>
         dplyr::group_by(dissimilarity_method, consistency.metric, celltype) |>
         dplyr::arrange(rate, .by_group = TRUE)
   } else if (group == "ident"){
      rsq <- df |> 
         dplyr::group_by(dissimilarity_method, consistency.metric, rate) |>
         dplyr::summarize(measure = mean(measure)) |>
         dplyr::group_by(dissimilarity_method, consistency.metric) |>
         dplyr::arrange(rate, .by_group = TRUE)
   } else {
      stop("Not supported grouping, it has to be either celltype or ident")
   }
   
   # remove NA or NaN
   rsq <- rsq[complete.cases(rsq),]
   
   type <- type[1] |> tolower()
   
   switch(type,
          "monotonic" = {
             rsq <- rsq |>
                dplyr::summarize(score_celltype = round(monotonicity_score(measure), 2))
             type.verb <- "Monotonic, expecting increasing score with rate.\n"
          },
          "referenceline" = {      
             rsq <- rsq <- rsq |>
                dplyr::summarize(score_celltype = round(fit.ReferenceLine(x = rate, y = measure)$r.squared, 2),
                                 pval = round(fit.ReferenceLine(x = rate, y = measure)$p.value, 2),
                                 p.val_fill = ifelse(pval < 0.05, "sig", "ns")
                )
             type.verb <- "fit.ReferenceLine, expecting a perfect curve of intercept = 0 and slope = 1\n"},
          "constant" = {
             rsq <- rsq <- rsq |>
                dplyr::summarize(score_celltype = round(fit.Constant(x = as.numeric(rate), y = measure)$`1-rss`, 2),
                                 pval = round(fit.Constant(x = as.numeric(rate), y = measure)$p.value, 2),
                                 p.val_fill = ifelse(pval < 0.05, "sig", "ns")
                )
             type.verb <- "fit.Constant, expecting a flat curve with slope = 0\n"},
          stop(type, "is not a supported scoring metric.")
   )
   
   rsq <- 
      rsq |>
      dplyr::group_by(dissimilarity_method, consistency.metric) |>
      dplyr::mutate(score_mean = round(mean(score_celltype, na.rm = TRUE), 2))
   
   
   if(verbose){message("\nEvaluating metric using ", type.verb)}
   
   if(!return.df){
      dts <- unique(df[["dissimilarity_method"]])
      pls <- list()
      
      range.x <- c(min(as.numeric(df[["rate"]])),
                   max(as.numeric(df[["rate"]]))
      )
      
      for(d in dts){
         dftmp <- df |>
            dplyr::filter(dissimilarity_method == d)
         rsqtmp <- rsq |>
            dplyr::filter(dissimilarity_method == d)
         
         ncol <- length(unique(dftmp[["consistency.metric"]]))
         
         plc <- dftmp |>
            ggplot2::ggplot(ggplot2::aes(rate, measure)) +
            ggplot2::geom_point(ggplot2::aes(color = celltype),
                                shape = 21,
                                alpha = 0.4) +
            ggplot2::geom_line(ggplot2::aes(color = celltype),
                               alpha = 0.4) +
            ggplot2::stat_summary(ggplot2::aes(group = consistency.metric),
                                  geom = "line",
                                  fun = function(y){mean(y, na.rm = TRUE)}) +
            ggplot2::stat_summary(ggplot2::aes(group = consistency.metric),
                                  geom = "point",
                                  fun = function(y){mean(y, na.rm = TRUE)}) +
            ggplot2::geom_label(data = rsqtmp,
                                ggplot2::aes(x = mean(range.x), y = 1.2,
                                             label = as.character(score_mean)),
                                alpha = 0.4,
                                show.legend = FALSE) +
            ggplot2::geom_hline(yintercept = 1,
                                linetype = "dotted") +
            ggplot2::facet_wrap(~consistency.metric,
                                ncol =  ncol,
                                drop = FALSE) +
            ggplot2::labs(x = xlabel,
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
