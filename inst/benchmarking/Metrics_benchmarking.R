
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



rand_shuffling_group <- function(vector,
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

rand_split_group <- function(vector,
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
                                  gene_list,
                                  reduction = TRUE,
                                  ndim = 30,
                                  normalization_method = "Log1p",
                                  dissimilarity_method = c("WasserStein", "Pseudobulk:Euclidean",
                                                           "Pseudobulk:Cosine", "Pseudobulk:Pearson",
                                                           "recip_classif:Match", "recip_classif:Score"),
                                  min_samples = 5,
                                  min_cells = 10,
                                  bparam = BiocParallel::SerialParam(),
                                  verbose = TRUE){
   
   scTypeEval <- run_processing_data(scTypeEval,
                                    ident = ident,
                                    sample = sample,
                                    normalization_method = normalization_method,
                                    min_samples = min_samples,
                                    min_cells = min_cells,
                                    verbose = verbose)
   scTypeEval <- add_gene_list(scTypeEval, gene_list)
   
   if(reduction){
      scTypeEval <- run_pca(scTypeEval,
                            ndim = ndim)
   }
   
   for(m in dissimilarity_method){
      if(verbose){message(">.  Running ", m, "\n")}
      scTypeEval <- run_dissimilarity(scTypeEval,
                                      reduction = reduction,
                                      method = m,
                                      bparam = bparam,
                                      verbose = verbose)
   }
   
   return(scTypeEval)
   
}

wrapper_plots <- function(scTypeEval,
                          reduction_slot = "all", 
                          dissimilarity_slot = "all",
                          reduction = TRUE,
                          label = TRUE,
                          dims = c(1,2),
                          show_legend = FALSE,
                          dir_path,
                          height =15,
                          width = 18,
                          ...){
   if(reduction){
      pca <- plot_pca(scTypeEval,
                      reduction_slot = reduction_slot,
                      label = label,
                      dims = dims,
                      show_legend = show_legend)
   } else {
      pca <- c()
   }
   
   mds <- plot_mds(scTypeEval,
                   dissimilarity_slot = dissimilarity_slot,
                   label = label,
                   dims = dims,
                   show_legend = show_legend)
   
   ph <- plot_heatmap(scTypeEval,
                      dissimilarity_slot = dissimilarity_slot,
                      ...)
   
   # Export as PDF
   dir <- paste0(dir_path, ".pdf")
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
wr_missclasify <- function(count_matrix,
                           metadata,
                           ident,
                           sample,
                           rates = c(1, 0.9, 0.75, 0.5, 0.25, 0), # proportion of cell labels kept
                           replicates = 3,
                           gene_list = NULL,
                           reduction = TRUE,
                           ndim = 30,
                           black_list = NULL,
                           dir = NULL,
                           normalization_method = "Log1p",
                           dissimilarity_method = c("WasserStein", "Pseudobulk:Euclidean",
                                                    "Pseudobulk:Cosine", "Pseudobulk:Pearson",
                                                    "recip_classif:Match", "recip_classif:Score"),
                           int_val_metric = c("silhouette", "NeighborhoodPurity",
                                             "ward_PropMatch", "Orbital_medoid",
                                             "Average_similarity", "2label_silhouette"),
                           min_samples = 5,
                           min_cells = 10,
                           ncores = 1,
                           bparam = NULL,
                           progressbar = FALSE,
                           seed = 22,
                           knn_graph_k = 5,
                           hclust_method = "ward.D2",
                           save_plots = TRUE,
                           verbose = TRUE){
   
   if(!is.null(dir)){
      if(verbose){message("\nResults will be stored at ", dir)}
      dir.create(dir, showWarnings = FALSE)
      
      if(save_plots){
         pca_dir <- file.path(dir, "Plots_Missclassify")
         dir.create(pca_dir, showWarnings = FALSE)
      }
   }
   
   # produce missclassification on metadata
   # get seed for each replicate
   sds <- seed + seq_len(replicates)
   names(sds) <- seq_len(replicates)
   
   original_vector <- metadata[[ident]]
   # vector for annotations
   annotations <- c()
   
   for(s in names(sds)){
      for(r in rates){
         ra <- 1-r # proportion of cells to shuffle
         new_vector <- rand_shuffling_group(vector = original_vector,
                                            group = metadata[[sample]],
                                            rate = ra,
                                            seed = sds[[s]])
         new_name <- paste("R", s, r, ident, sep = "_")
         metadata[[new_name]] <- new_vector
         annotations <- c(annotations, new_name)
         
      }
   }
   
   ## Create sc object
   sc <- create_scTypeEval(matrix = count_matrix,
                           metadata = metadata,
                           active_ident = ident,
                           black_list = black_list)
   # get the gene list, the same for every run
   if(is.null(gene_list)){
      sc_gl <- run_processing_data(sc,
                                  sample = sample,
                                  normalization_method = normalization_method,
                                  min_samples = min_samples,
                                  min_cells = min_cells,
                                  verbose = verbose)
      sc_gl <- run_hvg(sc_gl,
                       ncores = ncores)
      gl <- sc_gl@gene_lists
   } else {
      gl <- gene_list
   }
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   if(verbose){message("Running loop of annotations ")}
   df_res <- lapply(annotations,
                    function(ann){
                       tryCatch(
                          {
                             sc_tmp <- wrapper_dissimilarity(sc,
                                                             ident = ann,
                                                             sample = sample,
                                                             gene_list = gl,
                                                             reduction = reduction,
                                                             ndim = ndim,
                                                             normalization_method = normalization_method,
                                                             dissimilarity_method = dissimilarity_method,
                                                             min_samples = min_samples,
                                                             min_cells = min_cells,
                                                             bparam = param,
                                                             verbose = verbose
                             )
                             # data.frame with consistency outcome
                             res <- get_consistency(sc_tmp)
                             
                             # accommodate extra data
                             res <- res |>
                                dplyr::mutate(rate = as.numeric(as.character(strsplit(ann, "_")[[1]][3])),
                                              rep = strsplit(ann, "_")[[1]][2],
                                              original_ident = !!ident,
                                              task = "Missclassification"
                                )
                             
                             # render PCAs
                             # only produce for one replicate of the seeds
                             if(stringr::str_split(ann, "_")[[1]][2] == 1){
                                # save pdf if indicated
                                if(save_plots){
                                   if(verbose){message("\nProducing Plots for ", ann, "\n")}
                                   fp <- file.path(pca_dir, ann)
                                   wrapper_plots(sc_tmp,
                                                 dir_path = fp,
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
   df_res <- do.call(rbind, df_res)
   
   if(!is.null(dir)){
      saveRDS(df_res,
              file.path(dir,
                        paste0("Missclassification_", ident, ".rds")))
   }
   return(df_res)
   
}

# wrapper to evaluate the number of samples
wr_nsamples <- function(count_matrix,
                        metadata,
                        ident,
                        sample,
                        rates = c(1, 0.9, 0.7, 0.5), # proportion of cell labels kept
                        replicates = 3,
                        gene_list = NULL,
                        reduction = TRUE,
                        ndim = 30,
                        black_list = NULL,
                        dir = NULL,
                        normalization_method = "Log1p",
                        dissimilarity_method = c("WasserStein", "Pseudobulk:Euclidean",
                                                 "Pseudobulk:Cosine", "Pseudobulk:Pearson",
                                                 "recip_classif:Match", "recip_classif:Score"),
                        int_val_metric = c("silhouette", "NeighborhoodPurity",
                                          "ward_PropMatch", "Orbital_medoid",
                                          "Average_similarity", "2label_silhouette"),
                        min_samples = 5,
                        min_cells = 10,
                        ncores = 1,
                        bparam = NULL,
                        progressbar = FALSE,
                        seed = 22,
                        knn_graph_k = 5,
                        hclust_method = "ward.D2",
                        save_plots = TRUE,
                        verbose = TRUE){
   
   if(!is.null(dir)){
      if(verbose){message("\nResults will be stored at ", dir)}
      dir.create(dir, showWarnings = FALSE)
      
      if(save_plots){
         pca_dir <- file.path(dir, "Plots_Nsamples")
         dir.create(pca_dir, showWarnings = FALSE)
      }
   }
   
   
   # produce diferents combinations of samples
   # get seed for each replicate
   sds <- seed + seq_len(replicates)
   names(sds) <- seq_len(replicates)
   
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
   if(is.null(gene_list)){
      ## Create sc object
      sc <- create_scTypeEval(matrix = count_matrix,
                              metadata = metadata,
                              active_ident = ident,
                              black_list = black_list)
      sc_gl <- run_processing_data(sc,
                                  sample = sample,
                                  normalization_method = normalization_method,
                                  min_samples = min_samples,
                                  min_cells = min_cells,
                                  verbose = verbose)
      sc_gl <- run_hvg(sc_gl,
                       ncores = ncores)
                                         dplyr::mutate(rate = as.numeric(as.character(strsplit(ann, "_")[[1]][3])),
   } else {
      gl <- gene_list
   }
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   if(verbose){message("Running loop of annotations ")}
   df_res <- lapply(names(nss),
                    function(ns){
                       tryCatch(
                          {
                             ## run scTypeEval
                             # create scTypeEval object with the number of samples
                             md <- metadata[metadata[[sample]] %in% nss[[ns]],]
                             sc_tmp <- create_scTypeEval(matrix = count_matrix[,rownames(md)],
                                                         metadata = md,
                                                         gene_lists = gl,
                                                         black_list = black_list)
                             
                             sc_tmp <- wrapper_dissimilarity(sc_tmp,
                                                             ident = ident,
                                                             sample = sample,
                                                             gene_list = gl,
                                                             reduction = reduction,
                                                             ndim = ndim,
                                                             normalization_method = normalization_method,
            df_res <- do.call(rbind, df_res)
                                                             min_samples = min_samples,
                                                             min_cells = min_cells,
                                                             bparam = param,
                                                             verbose = verbose
                             )
                             # data.frame with consistency outcome
                             res <- get_consistency(sc_tmp)
                             
                             # accommodate extra data
                             res <- res |>
                                dplyr::mutate(rate = as.numeric(as.character(strsplit(ns, "_")[[1]][3])),
                                              rep = strsplit(ns, "_")[[1]][2],
                                              original_ident = !!ident,
                                              task = "NSamples"
                                )
                             
                             # render PCAs
                             # only produce for one replicate of the seeds
                             if(stringr::str_split(ns, "_")[[1]][2] == 1){
                                # save pdf if indicated
                                if(save_plots){
                                   if(verbose){message("\nProducing Plots for ", ns, "\n")}
                                   fp <- file.path(pca_dir, ns)
                                   wrapper_plots(sc_tmp,
                                                 dir_path = fp,
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
   df_res <- do.call(rbind, df_res)
   
   if(!is.null(dir)){
      saveRDS(df_res,
              file.path(dir,
                        paste0("NSamples_", ident, ".rds")))
   }
   
   return(df_res)
   
}

# wrapper to evaluate the consistency metrics when excluding some cell types
wr_nct <- function(count_matrix,
                   metadata,
                   ident,
                   sample,
                   gene_list = NULL,
                   reduction = TRUE,
                   ndim = 30,
                   black_list = NULL,
                   dir = NULL,
                   normalization_method = "Log1p",
                   dissimilarity_method = c("WasserStein", "Pseudobulk:Euclidean",
                                            "Pseudobulk:Cosine", "Pseudobulk:Pearson",
                                            "recip_classif:Match", "recip_classif:Score"),
                   int_val_metric = c("silhouette", "NeighborhoodPurity",
                                     "ward_PropMatch", "Orbital_medoid",
                                     "Average_similarity", "2label_silhouette"),
                   min_samples = 5,
                   min_cells = 10,
                   ncores = 1,
                   bparam = NULL,
                   progressbar = FALSE,
                   seed = 22,
                   knn_graph_k = 5,
                   hclust_method = "ward.D2",
                   save_plots = TRUE,
                   verbose = TRUE,
                   down_sample_comb = 30){
   
   if(!is.null(dir)){
      if(verbose){message("Results will be stored at ", dir)}
      dir.create(dir, showWarnings = FALSE)
      
      if(save_plots){
         pca_dir <- file.path(dir, "Plots_Nct")
         dir.create(pca_dir, showWarnings = FALSE)
      }
   }
   
   
   # get all possible combinations
   allcts <- unique(metadata[[ident]])
   allcts <- allcts[!is.na(allcts)]
   
   
   # Generate combinations of lengths 2 to n-1
   cts <- sample_variable_length_combinations(allcts,
                                              num_samples = down_sample_comb,
                                              seed = seed)
   # add all cell types
   allelements <- paste(allcts, collapse = "-")
   cts[[allelements]] <- allcts
   
   # get the gene list, the same for every run
   if(is.null(gene_list)){
      ## Create sc object
      sc <- create_scTypeEval(matrix = count_matrix,
                              metadata = metadata,
                              active_ident = ident,
                              black_list = black_list)
      sc_gl <- run_processing_data(sc,
                                  sample = sample,
                                  normalization_method = normalization_method,
                                  min_samples = min_samples,
                                  min_cells = min_cells,
                                  verbose = verbose)
      sc_gl <- run_hvg(sc_gl,
                       ncores = ncores)
      gl <- sc_gl@gene_lists
   } else {
      gl <- gene_list
   }
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   df_res <- lapply(names(cts),
                    function(ns){
                       tryCatch(
                          {
                             # create scTypeEval object with the number of samples
                             md <- metadata[metadata[[ident]] %in% cts[[ns]],]
                             sc_tmp <- create_scTypeEval(matrix = count_matrix[,rownames(md)],
                                                         metadata = md,
                                                         gene_lists = gl,
                                                         black_list = black_list)
                             sc_tmp <- wrapper_dissimilarity(sc_tmp,
                                                             ident = ident,
                                                             sample = sample,
                                                             gene_list = gl,
                                                             reduction = reduction,
                                                             ndim = ndim,
                                                             normalization_method = normalization_method,
                                                             dissimilarity_method = dissimilarity_method,
                                                             min_samples = min_samples,
                                                             min_cells = min_cells,
                                                             bparam = param,
                                                             verbose = verbose
                             )
                             # data.frame with consistency outcome
                             res <- get_consistency(sc_tmp)
                             
                             # accommodate extra data
                             res <- res |>
                                dplyr::mutate(rate = ns,
                                              rep = NA,
                                              original_ident = !!ident,
                                              task = "Nct"
                                )
                             
                             # render plots
                             # save pdf if indicated
                             if(save_plots){
                                if(verbose){message("\nProducing Plots for ", ns, "\n")}
                                fp <- file.path(pca_dir, ns)
                                wrapper_plots(sc_tmp,
                                              dir_path = fp,
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
   df_res <- do.call(rbind, df_res) |>
      # convert to factor to then compute fit_constant
      dplyr::mutate(rate = factor(rate,
                                  levels = names(cts))
      )
   
   if(!is.null(dir)){
      saveRDS(df_res,
              file.path(dir,
                        paste0("Nct_", ident, ".rds")))
   }
   
   return(df_res)
}

# wrapper to evaluate the number of cells in a cell type
wr_ncell <- function(count_matrix,
                     metadata,
                     ident,
                     sample,
                     rates = c(1, 0.75, 0.5, 0.25), # proportion of total samples kept
                     ctype = NULL, # which cell type lower the proportion
                     replicates = 3,
                     gene_list = NULL,
                     reduction = TRUE,
                     ndim = 30,
                     black_list = NULL,
                     dir = NULL,
                     normalization_method = "Log1p",
                     dissimilarity_method = c("WasserStein", "Pseudobulk:Euclidean",
                                              "Pseudobulk:Cosine", "Pseudobulk:Pearson",
                                              "recip_classif:Match", "recip_classif:Score"),
                     int_val_metric = c("silhouette", "NeighborhoodPurity",
                                       "ward_PropMatch", "Orbital_medoid",
                                       "Average_similarity", "2label_silhouette"),
                     min_samples = 5,
                     min_cells = 10,
                     ncores = 1,
                     bparam = NULL,
                     progressbar = FALSE,
                     seed = 22,
                     knn_graph_k = 5,
                     hclust_method = "ward.D2",
                     save_plots = TRUE,
                     verbose = TRUE
){
   
   if(!is.null(dir)){
      if(verbose){message("\nResults will be stored at ", dir)}
      dir.create(dir, showWarnings = FALSE)
      
      if(save_plots){
         pca_dir <- file.path(dir, "Plots_NCell")
         dir.create(pca_dir, showWarnings = FALSE)
      }
   }
   
   original_vector <- as.character(metadata[[ident]])
   
   # get cell type to split, defined or most abundant
   if(is.null(ctype)){
      ctype <- names(which.max(table(original_vector)))
   }
   if(verbose){message("\nReducing number of cells for ", ctype, "\n")}
   
   
   # produce diferents combinations of samples
   # get seed for each replicate
   sds <- seed + seq_len(replicates)
   names(sds) <- seq_len(replicates)
   
   # get number of cells from ctype
   nb <- table(metadata[[ident]])
   nb <- nb[ctype]
   
   
   # get the gene list, the same for every run
   if(is.null(gene_list)){
      ## Create sc object
      sc <- create_scTypeEval(matrix = count_matrix,
                              metadata = metadata,
                              active_ident = ident,
                              black_list = black_list)
      sc_gl <- run_processing_data(sc,
                                  sample = sample,
                                  normalization_method = normalization_method,
                                  min_samples = min_samples,
                                  min_cells = min_cells,
                                  verbose = verbose)
      sc_gl <- run_hvg(sc_gl,
                       ncores = ncores)
      gl <- sc_gl@gene_lists
   } else {
      gl <- gene_list
   }
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   combi <- expand.grid(names(sds), rates)
   
   df_res <- lapply(seq_len(nrow(combi)),
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
                             sc_tmp <- create_scTypeEval(matrix = count_matrix[,rownames(md)],
                                                         metadata = md,
                                                         gene_lists = gl,
                                                         black_list = black_list)
                             
                             sc_tmp <- wrapper_dissimilarity(sc_tmp,
                                                             ident = ident,
                                                             sample = sample,
                                                             gene_list = gl,
                                                             reduction = reduction,
                                                             ndim = ndim,
                                                             normalization_method = normalization_method,
                                                             dissimilarity_method = dissimilarity_method,
                                                             min_samples = min_samples,
                                                             min_cells = min_cells,
                                                             bparam = param,
                                                             verbose = verbose
                             )
                             # data.frame with consistency outcome
                             res <- get_consistency(sc_tmp)
                             
                             # accommodate extra data
                             res <- res |>
                                dplyr::mutate(rate = as.numeric(as.character(ns)),
                                              rep = s,
                                              original_ident = !!ident,
                                              task = "NCell"
                                )
                             
                             # render PCAs
                             # only produce for one replicate of the seeds
                             if(s == 1){
                                # save pdf if indicated
                                # save pdf if indicated
                                if(save_plots){
                                   if(verbose){message("\nProducing Plots for ", ns, "\n")}
                                   fp <- file.path(pca_dir, ns)
                                   wrapper_plots(sc_tmp,
                                                 dir_path = fp,
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
   df_res <- do.call(rbind, df_res)
   ctype_name <- purge_label(ctype)
   
   if(!is.null(dir)){
      saveRDS(df_res,
              file.path(dir,
                        paste0("NCell_", ctype_name, "_", ident, ".rds")))
   }
   
   
   return(df_res)
   
}



# wrapper to hirarchally merge cell types
wr_merge_ct <- function(count_matrix,
                       metadata,
                       ident,
                       sample,
                       distance_method = "euclidean",
                       gene_list = NULL,
                       reduction = TRUE,
                       ndim = 30,
                       black_list = NULL,
                       dir = NULL,
                       normalization_method = "Log1p",
                       dissimilarity_method = c("WasserStein", "Pseudobulk:Euclidean",
                                                "Pseudobulk:Cosine", "Pseudobulk:Pearson",
                                                "recip_classif:Match", "recip_classif:Score"),
                       int_val_metric = c("silhouette", "NeighborhoodPurity",
                                         "ward_PropMatch", "Orbital_medoid",
                                         "Average_similarity", "2label_silhouette"),
                       min_samples = 5,
                       min_cells = 10,
                       ncores = 1,
                       bparam = NULL,
                       progressbar = FALSE,
                       knn_graph_k = 5,
                       hclust_method = "ward.D2",
                       save_plots = TRUE,
                       verbose = TRUE){
   
   if(!is.null(dir)){
      if(verbose){message("\nResults will be stored at ", dir)}
      dir.create(dir, showWarnings = FALSE)
      
      if(save_plots){
         pca_dir <- file.path(dir, "Plots_mergeCT")
         dir.create(pca_dir, showWarnings = FALSE)
      }
   }
   
   ## Create sc object
   sc <- create_scTypeEval(matrix = count_matrix,
                           metadata = metadata,
                           active_ident = ident,
                           black_list = black_list)
   
   sc <- run_processing_data(sc,
                            sample = sample,
                            normalization_method = normalization_method,
                            min_samples = min_samples,
                            min_cells = min_cells,
                            verbose = verbose)
   if(is.null(gene_list)){
      
      sc <- run_hvg(sc,
                    ncores = ncores)
      gl <- sc@gene_lists
   } else {
      gl <- gene_list
   }
   
   gene_list <- check_genelist(sc, gene_list, verbose = verbose)
   black_list <- check_blacklist(sc, black_list, verbose = verbose)
   mat_ident <- general_filtering(sc@data[["single-cell"]],
                                   black_list = black_list,
                                   gene_list = gene_list,
                                   verbose = verbose)
   mat <- mat_ident@matrix
   
   idents <- as.character(mat_ident@ident)
   # keep only celltypes fulfulling min_samples and min_cells conditions
   metadata <- metadata[colnames(mat),]
   
   original_vector <- factor(purge_label(metadata[[ident]]))
   j <- length(unique(original_vector))
   new_name <- paste("J", j, ident, sep = "_")
   metadata[[new_name]] <- original_vector
   annotations <- new_name
   
   while(length(unique(idents))>2){
      centroids <- compute_centroids(norm_mat = mat,
                                     ident = idents)
      dist <- get_distance(norm_mat = centroids,
                           distance_method = distance_method,
                           verbose = FALSE)
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
   sc <- create_scTypeEval(matrix = count_matrix[,rownames(metadata)],
                           metadata = metadata,
                           gene_lists = gl,
                           black_list = black_list)
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   df_res <- lapply(annotations,
                    function(ann){
                       tryCatch(
                          {
                             sc_tmp <- wrapper_dissimilarity(sc,
                                                             ident = ann,
                                                             sample = sample,
                                                             gene_list = gl,
                                                             reduction = reduction,
                                                             ndim = ndim,
                                                             normalization_method = normalization_method,
                                                             dissimilarity_method = dissimilarity_method,
                                                             min_samples = min_samples,
                                                             min_cells = min_cells,
                                                             bparam = param,
                                                             verbose = verbose
                             )
                             # data.frame with consistency outcome
                             res <- get_consistency(sc_tmp)
                             
                             # accommodate extra data
                             res <- res |>
                                dplyr::mutate(
                                   rate = as.numeric(as.character(strsplit(ann, "_")[[1]][2])),
                                   rep = NA,
                                   original_ident = !!ident,
                                   task = "mergeCT"
                                )
                             
                             # render PCAs
                             # save pdf if indicated
                             if(save_plots){
                                if(verbose){message("\nProducing Plots for ", ann, "\n")}
                                fp <- file.path(pca_dir, ann)
                                wrapper_plots(sc_tmp,
                                              dir_path = fp,
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
   df_res <- do.call(rbind, df_res)
   
   if(!is.null(dir)){
      saveRDS(df_res,
              file.path(dir,
                        paste0("mergeCT_", ident, ".rds")))
   }
   
   return(df_res)
   
}

# wrapper to split a cell type into 2 different, but highly similar transcriptomically
wr_split_cell_type <- function(count_matrix,
                             metadata,
                             ident,
                             sample,
                             rates = c(1, 0.5),
                             ctype = NULL, # which cell type to split
                             replicates = 3,
                             gene_list = NULL,
                             reduction = TRUE,
                             ndim = 30,
                             black_list = NULL,
                             dir = NULL,
                             normalization_method = "Log1p",
                             dissimilarity_method = c("WasserStein", "Pseudobulk:Euclidean",
                                                      "Pseudobulk:Cosine", "Pseudobulk:Pearson",
                                                      "recip_classif:Match", "recip_classif:Score"),
                             int_val_metric = c("silhouette", "NeighborhoodPurity",
                                               "ward_PropMatch", "Orbital_medoid",
                                               "Average_similarity", "2label_silhouette"),
                             min_samples = 5,
                             min_cells = 10,
                             ncores = 1,
                             bparam = NULL,
                             progressbar = FALSE,
                             seed = 22,
                             knn_graph_k = 5,
                             hclust_method = "ward.D2",
                             save_plots = TRUE,
                             verbose = TRUE){
   
   if(!is.null(dir)){
      if(verbose){message("\nResults will be stored at ", dir)}
      dir.create(dir, showWarnings = FALSE)
      
      if(save_plots){
         pca_dir <- file.path(dir, "Plots_SplitCellType")
         dir.create(pca_dir, showWarnings = FALSE)
      }
   }
   
   # produce missclassification on metadata
   # get seed for each replicate
   sds <- seed + seq_len(replicates)
   names(sds) <- seq_len(replicates)
   
   original_vector <- as.character(metadata[[ident]])
   
   # get cell type to split, defined or most abundant
   if(is.null(ctype)){
      ctype <- names(which.max(table(original_vector)))
   }
   if(verbose){message("\nReducing number of cells for ", ctype, "\n")}
   
   
   # vector for annotations
   annotations <- c()
   
   for(s in names(sds)){
      for(r in rates){
         ra <- 1-r # proportion of cells to shuffle
         new_vector <- rand_split_group(vector = original_vector,
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
   sc <- create_scTypeEval(matrix = count_matrix,
                           metadata = metadata,
                           active_ident = ident,
                           black_list = black_list)
   # get the gene list, the same for every run
   if(is.null(gene_list)){
      sc_gl <- run_processing_data(sc,
                                  sample = sample,
                                  normalization_method = normalization_method,
                                  min_samples = min_samples,
                                  min_cells = min_cells,
                                  verbose = verbose)
      sc_gl <- run_hvg(sc_gl,
                       ncores = ncores)
      gl <- sc_gl@gene_lists
   } else {
      gl <- gene_list
   }
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   if(verbose){message("Running loop of annotations ")}
   df_res <- lapply(annotations,
                    function(ann){
                       tryCatch(
                          {
                             sc_tmp <- wrapper_dissimilarity(sc,
                                                             ident = ann,
                                                             sample = sample,
                                                             gene_list = gl,
                                                             reduction = reduction,
                                                             ndim = ndim,
                                                             normalization_method = normalization_method,
                                                             dissimilarity_method = dissimilarity_method,
                                                             min_samples = min_samples,
                                                             min_cells = min_cells,
                                                             bparam = param,
                                                             verbose = verbose
                             )
                             # data.frame with consistency outcome
                             res <- get_consistency(sc_tmp)
                             
                             # accommodate extra data
                             res <- res |>
                                dplyr::mutate(rate = as.numeric(as.character(strsplit(ann, "_")[[1]][3])),
                                              rep = strsplit(ann, "_")[[1]][2],
                                              original_ident = !!ident,
                                              task = "SplitCelltype"
                                )
                             
                             # render PCAs
                             # only produce for one replicate of the seeds
                             if(stringr::str_split(ann, "_")[[1]][2] == 1){
                                # save pdf if indicated
                                if(save_plots){
                                   if(verbose){message("\nProducing Plots for ", ann, "\n")}
                                   fp <- file.path(pca_dir, ann)
                                   wrapper_plots(sc_tmp,
                                                 dir_path = fp,
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
   df_res <- do.call(rbind, df_res)
   ctype_name <- purge_label(ctype)
   
   if(!is.null(dir)){
      saveRDS(df_res,
              file.path(dir,
                        paste0("SplitCelltype_", ctype_name, "_", ident, ".rds")))
   }
   
   return(df_res)
}


wr_assay_plot <- function(df,
                         type = c("Monotonic", "ReferenceLine", "Constant", "drop"),
                         group = c("celltype", "ident"),
                         return_df = FALSE,
                         verbose = TRUE,
                         xlabel = "rate",
                         title = "",
                         combine = TRUE){
   
   group <- group[1] |> tolower()
   
   if(group == "celltype"){
      rsq <- df |> 
         dplyr::group_by(dissimilarity_method, consistency_metric, celltype, rate) |>
         dplyr::summarize(measure = mean(measure)) |>
         dplyr::group_by(dissimilarity_method, consistency_metric, celltype) |>
         dplyr::arrange(rate, .by_group = TRUE)
   } else if (group == "ident"){
      rsq <- df |> 
         dplyr::group_by(dissimilarity_method, consistency_metric, rate) |>
         dplyr::summarize(measure = mean(measure)) |>
         dplyr::group_by(dissimilarity_method, consistency_metric) |>
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
             type_verb <- "Monotonic, expecting increasing score with rate.\n"
          },
          "referenceline" = {      
             rsq <- rsq <- rsq |>
                dplyr::summarize(score_celltype = round(fit_reference_line(x = rate, y = measure)$r_squared, 2),
                                 pval = round(fit_reference_line(x = rate, y = measure)$p_value, 2),
                                 p_val_fill = ifelse(pval < 0.05, "sig", "ns")
                )
             type_verb <- "fit_reference_line, expecting a perfect curve of intercept = 0 and slope = 1\n"},
          "constant" = {
             rsq <- rsq <- rsq |>
                dplyr::summarize(score_celltype = round(fit_constant(x = as.numeric(rate), y = measure)$one_minus_rss, 2),
                                 pval = round(fit_constant(x = as.numeric(rate), y = measure)$p_value, 2),
                                 p_val_fill = ifelse(pval < 0.05, "sig", "ns")
                )
             type_verb <- "fit_constant, expecting a flat curve with slope = 0\n"},
          "drop" = {
             rsq <- rsq <- rsq |>
                dplyr::summarize(score_celltype = round(consistency_drop(x = measure, y = rate), 2))
             type_verb <- "fit Consistency drop, expecting a drop in consistency in smaller rates.\n"
          },
          stop(type, " is not a supported scoring metric.")
   )
   
   rsq <- 
      rsq |>
      dplyr::group_by(dissimilarity_method, consistency_metric) |>
      dplyr::mutate(score_mean = round(mean(score_celltype, na.rm = TRUE), 2))
   
   
   if(verbose){message("\nEvaluating metric using ", type_verb)}
   
   if(!return_df){
      dts <- unique(df[["dissimilarity_method"]])
      pls <- list()
      
      range_x <- c(min(as.numeric(df[["rate"]])),
                   max(as.numeric(df[["rate"]]))
      )
      
      for(d in dts){
         dftmp <- df |>
            dplyr::filter(dissimilarity_method == d)
         rsqtmp <- rsq |>
            dplyr::filter(dissimilarity_method == d)
         
         ncol <- length(unique(dftmp[["consistency_metric"]]))
         
         plc <- dftmp |>
            ggplot2::ggplot(ggplot2::aes(rate, measure)) +
            ggplot2::geom_point(ggplot2::aes(color = celltype),
                                shape = 21,
                                alpha = 0.4) +
            ggplot2::geom_line(ggplot2::aes(color = celltype),
                               alpha = 0.4) +
            ggplot2::stat_summary(ggplot2::aes(group = consistency_metric),
                                  geom = "line",
                                  fun = function(y){mean(y, na.rm = TRUE)}) +
            ggplot2::stat_summary(ggplot2::aes(group = consistency_metric),
                                  geom = "point",
                                  fun = function(y){mean(y, na.rm = TRUE)}) +
            ggplot2::geom_label(data = rsqtmp,
                                ggplot2::aes(x = mean(range_x), y = 1.2,
                                             label = as.character(score_mean)),
                                alpha = 0.4,
                                show_legend = FALSE) +
            ggplot2::geom_hline(yintercept = 1,
                                linetype = "dotted") +
            ggplot2::facet_wrap(~consistency_metric,
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
