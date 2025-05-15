downsample_factors <- function(df, factor_col, threshold, seed = 22) {
   set.seed(seed) # Set seed for reproducibility
   
   # Group by the factor column and apply downsampling
   df_downsampled <- df |> 
      tibble::rownames_to_column("cellid") |> 
      dplyr::group_by(!!rlang::sym(factor_col)) |> 
      dplyr::group_map(~ { 
         if (nrow(.x) > threshold) { 
            .x[sample(nrow(.x), threshold), ] # Randomly sample rows if exceeding threshold 
         } else { 
            .x # Keep all rows if below threshold 
         } 
      }, .keep = TRUE) |> # Ensure the grouping column is kept 
      dplyr::bind_rows() |> 
      dplyr::ungroup() |> 
      tibble::column_to_rownames("cellid")
   
   return(df_downsampled)
}

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

rand.shuffling <- function(vector,
                           rate = 0.1, # proportion of labels to shuffle
                           prob = NULL,
                           seed = 22) {
   set.seed(seed)
   
   # Calculate the number of labels to shuffle
   n <- length(vector)
   shuffle_n <- round(n * rate)
   
   # Get the proportions of cells
   proportions_cells <- table(vector) |> prop.table()
   
   # Indices to shuffle
   shuffle_indices <- sample(n, shuffle_n)
   
   # Extract the values to shuffle
   original_values <- vector[shuffle_indices]
   
   # Ensure no value remains in its original label
   shuffled_values <- original_values
   while (any(shuffled_values == original_values)) {
      ind_diff <- shuffled_values == original_values
      if (is.null(prob)) {
         shuffled_values[ind_diff] <- sample(vector,
                                             length(shuffled_values[ind_diff]),
                                             replace = TRUE)
      } else {
         # Create a matrix of probabilities for each value
         prob_matrix <- sapply(unique(vector), function(l) {
            prob.custom <- prob[, as.character(l)]
            prob.custom <- 1 / prob.custom
            prob.custom[is.infinite(prob.custom)] <- 0
            # Get weights based on distance to centroids and proportions
            prob.custom <- prob.custom * proportions_cells
            # Normalize probabilities
            prob.custom <- prob.custom / sum(prob.custom)
            prob.custom
         })
         
         colnames(prob_matrix) <- unique(vector)
         
         svals <- shuffled_values[ind_diff]
         
         for (l in unique(svals)) {
            svals[svals == l] <- sample(unique(vector),
                                        length(svals[svals == l]),
                                        replace = TRUE,
                                        prob = prob_matrix[, l])
         }
         
         shuffled_values[ind_diff] <- svals
      }
   }
   
   # Replace original values with shuffled values
   vector[shuffle_indices] <- shuffled_values
   
   return(vector)
}

rand.shuffling2 <- function(vector,
                            rate = 0.1, # proportion of labels to shuffle
                            seed = 22) {
   set.seed(seed)
   
   nvector <- vector
   # Calculate the number of labels to shuffle
   n <- length(vector)
   shuffle_n <- floor(round(n * rate))
   its <- sample(1:n, shuffle_n)  # instances to shuffle
   
   nvector[its] <- sample(vector, shuffle_n)
   
   return(nvector)
}

df.missclassify <- function(metadata,
                            annotations,
                            rates,
                            prefix = "R-") {
   for (a in annotations) {
      original_vector <- metadata[[a]]
      
      for (r in rates) {
         ra <- 1 - r # proportion of cells to shuffle
         new_vector <- rand.shuffling(original_vector,
                                      rate = ra)
         new_name <- paste0("R-", r, "_", a)
         metadata[[new_name]] <- new_vector
         
      }
   }
   
   return(metadata)
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


# wrapper to evaluate the missclassifiction rate 
wr.missclasify <- function(count_matrix,
                           metadata,
                           ident,
                           sample,
                           rates = c(1, 0.9, 0.75, 0.5, 0.25, 0), # proportion of cell labels kept
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
                           save.PCA = TRUE,
                           dims = c(1,2),
                           show.legend = FALSE,
                           label = TRUE,
                           verbose = TRUE){
   
   if(!is.null(dir)){
      if(verbose){message("\nResults will be stored at ", dir)}
      dir.create(dir)
      
      if(save.PCA){
         pca.dir <- file.path(dir, "PCA_Missclassify")
         dir.create(pca.dir)
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
         new_vector <- rand.shuffling2(original_vector,
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
   
   df.res <- list()
   pcas <- list()
   ## run scTypeEval
   for(ann in annotations){
      res <- Run.scTypeEval(scTypeEval = sc,
                            ident = ann,
                            ident_GroundTruth = ident,
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
                            ncores = ncores,
                            bparam = bparam,
                            progressbar = progressbar,
                            verbose = verbose)
      
      # accommodate extra data
      res <- res |>
         dplyr::mutate(rate = as.numeric(as.character(strsplit(ns, "_")[[1]][3])),
                       rep = strsplit(ann, "_")[[1]][2],
                       original.ident = ident,
                       task = "Missclassification"
         )
      df.res[[ann]] <- res
      
      # render PCAs
      # only produce for one replicate of the seeds
      if(stringr::str_split(ann, "_")[[1]][2] != 1) {next}
      # save pdf if indicated
      if(save.PCA){
         
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
                       ncores = ncores,
                       bparam = bparam,
                       progressbar = progressbar,
                       verbose = verbose)
         pcas[[ann]] <- pcs
         
         for (p in names(pcs)){
            pdf(file.path(pca.dir, paste0(p, "_", ann, ".pdf")),
                width = 4, height = 4)
            print(pcs[[p]])
            dev.off()
         }
         
      }
      
   }
   
   # concatenate all results
   df.res <- do.call(rbind, df.res)
   
   if(!is.null(dir)){
      saveRDS(df.res,
              file.path(dir,
                        paste0("Missclassification_", ident, ".rds")))
   }
   
   ret.list <- list(df.res = df.res,
                    PCA = pcas)
   
   return(ret.list)
   
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
                        save.PCA = TRUE,
                        dims = c(1,2),
                        show.legend = FALSE,
                        label = TRUE,
                        verbose = TRUE){
   
   if(!is.null(dir)){
      if(verbose){message("\nResults will be stored at ", dir)}
      dir.create(dir)
      
      if(save.PCA){
         pca.dir <- file.path(dir, "PCA_NSamples")
         dir.create(pca.dir)
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
   
   df.res <- list()
   pcas <- list()
   ## run scTypeEval
   for(ns in names(nss)){
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
                            ncores = ncores,
                            bparam = bparam,
                            progressbar = progressbar,
                            verbose = verbose)
      
      # accommodate extra data
      res <- res |>
         dplyr::mutate(rate = as.numeric(as.character(strsplit(ns, "_")[[1]][3])),
                       rep = strsplit(ns, "_")[[1]][2],
                       original.ident = ident,
                       task = "NSamples"
         )
      df.res[[ns]] <- res
      
      # render PCAs
      # only produce for one replicate of the seeds
      if(strsplit(ns, "_")[[1]][2] != 1) {next}
      # save pdf if indicated
      if(save.PCA){
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
                       ncores = ncores,
                       bparam = bparam,
                       progressbar = progressbar,
                       verbose = verbose)
         pcas[[ns]] <- pcs
         
         
         for (p in names(pcs)){
            pdf(file.path(pca.dir, paste0(p, "_", ns, ".pdf")),
                width = 4, height = 4)
            print(pcs[[p]])
            dev.off()
         }
         
      }
      
   }
   
   # concatenate all results
   df.res <- do.call(rbind, df.res)
   
   if(!is.null(dir)){
      saveRDS(df.res,
              file.path(dir,
                        paste0("NSamples_", ident, ".rds")))
   }
   
   ret.list <- list(df.res = df.res,
                    PCA = pcas)
   
   return(ret.list)
   
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
                   save.PCA = TRUE,
                   dims = c(1,2),
                   show.legend = FALSE,
                   label = TRUE,
                   verbose = TRUE){
   
   if(!is.null(dir)){
      if(verbose){message("Results will be stored at ", dir)}
      dir.create(dir)
      
      if(save.PCA){
         pca.dir <- file.path(dir, "PCA_Nct")
         dir.create(pca.dir)
      }
   }
   
   
   # get all possible combinations
   allcts <- unique(metadata[[ident]])
   
   # Generate combinations of lengths 2 to n-1
   cts <- lapply(2:length(allcts),
                 function(k) combn(allcts, k, simplify = FALSE))
   # Flatten the list into a single vector of combinations
   cts <- do.call(c, cts)
   
   # name list
   names(cts) <- lapply(cts, paste, collapse = "-")
   
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
   
   df.res <- list()
   pcas <- list()
   ## run scTypeEval
   for(ns in names(cts)){
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
                                  ncores = ncores,
                                  bparam = bparam,
                                  progressbar = progressbar,
                                  verbose = verbose)
            
            # accommodate extra data
            res <- res |>
               dplyr::mutate(rate = ns,
                             rep = NA,
                             original.ident = ident,
                             task = "NSamples"
               )
            df.res[[ns]] <- res
            
            # render PCAs
            # only produce for one replicate of the seeds
            # save pdf if indicated
            if(save.PCA){
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
                             ncores = ncores,
                             bparam = bparam,
                             progressbar = progressbar,
                             verbose = verbose)
               pcas[[ns]] <- pcs
               
               
               for (p in names(pcs)){
                  pdf(file.path(pca.dir, paste0(p, "_", ns, ".pdf")),
                      width = 4, height = 4)
                  print(pcs[[p]])
                  dev.off()
               }
               
            }
         },
         error = function(e){
            message("\n-X- scTypeEval failed for", ns, "\n")
            message(e)
         }
      )
      
   }
   
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
   
   ret.list <- list(df.res = df.res,
                    PCA = pcas)
   
   return(ret.list)
   
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
                     save.PCA = TRUE,
                     dims = c(1,2),
                     show.legend = FALSE,
                     label = TRUE,
                     verbose = TRUE){
   
   if(!is.null(dir)){
      if(verbose){message("\nResults will be stored at ", dir)}
      dir.create(dir)
      
      if(save.PCA){
         pca.dir <- file.path(dir, "PCA_NCell")
         dir.create(pca.dir)
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
   
   df.res <- list()
   pcas <- list()
   ## run scTypeEval
   for(s in names(sds)){
      for(ns in rates){
         # create scTypeEval object with the number of samples
         # subset cell type
         nn <- paste(s, ns, ident, sep = "_")
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
                               ncores = ncores,
                               bparam = bparam,
                               progressbar = progressbar,
                               verbose = verbose)
         
         # accommodate extra data
         res <- res |>
            dplyr::mutate(rate = as.numeric(as.character(ns)),
                          rep = s,
                          original.ident = ident,
                          task = "NCell"
            )
         df.res[[nn]] <- res
         
         # render PCAs
         # only produce for one replicate of the seeds
         if(s != 1) {next}
         # save pdf if indicated
         if(save.PCA){
            if(verbose){message("\nProducing PCAs for ", nn, "\n")}
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
                          ncores = ncores,
                          bparam = bparam,
                          progressbar = progressbar,
                          verbose = verbose)
            pcas[[nn]] <- pcs
            
            
            for (p in names(pcs)){
               pdf(file.path(pca.dir, paste0(p, "_", nn, ".pdf")),
                   width = 4, height = 4)
               print(pcs[[p]])
               dev.off()
            }
            
         }
         
      }
   }
   
   # concatenate all results
   df.res <- do.call(rbind, df.res)
   
   if(!is.null(dir)){
      saveRDS(df.res,
              file.path(dir,
                        paste0("NCell_", ctype, "_", ident, ".rds")))
   }
   
   ret.list <- list(df.res = df.res,
                    PCA = pcas)
   
   return(ret.list)
   
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
                       save.PCA = TRUE,
                       dims = c(1,2),
                       show.legend = FALSE,
                       label = TRUE,
                       verbose = TRUE){
   
   if(!is.null(dir)){
      if(verbose){message("\nResults will be stored at ", dir)}
      dir.create(dir)
      
      if(save.PCA){
         pca.dir <- file.path(dir, "PCA_mergeCT")
         dir.create(pca.dir)
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
   
   df.res <- list()
   pcas <- list()
   ## run scTypeEval
   for(ann in annotations){
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
                            ncores = ncores,
                            bparam = bparam,
                            progressbar = progressbar,
                            verbose = verbose)
      
      # accommodate extra data
      res <- res |>
         dplyr::mutate(
            rate = as.numeric(as.character(strsplit(ns, "_")[[1]][2])),
            rep = NA,
            original.ident = ident,
            task = "mergeCT"
         )
      df.res[[ann]] <- res
      
      # render PCAs
      # save pdf if indicated
      if(save.PCA){
         
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
                       ncores = ncores,
                       bparam = bparam,
                       progressbar = progressbar,
                       verbose = verbose)
         pcas[[ann]] <- pcs
         
         
         for (p in names(pcs)){
            pdf(file.path(pca.dir, paste0(p, "_", ann, ".pdf")),
                width = 4, height = 4)
            print(pcs[[p]])
            dev.off()
         }
         
      }
      
   }
   
   # concatenate all results
   df.res <- do.call(rbind, df.res)
   
   if(!is.null(dir)){
      saveRDS(df.res,
              file.path(dir,
                        paste0("mergeCT_", ident, ".rds")))
   }
   
   ret.list <- list(df.res = df.res,
                    PCA = pcas)
   
   return(ret.list)
   
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
                                 trim = trim)) |>
      dplyr::group_by(gene.list, data.type, consistency.metric)
   
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
            ggplot2::geom_point(ggplot2::aes(color = celltype),
                                shape = 1) +
            ggplot2::stat_summary(ggplot2::aes(group = consistency.metric),
                                  geom = "line",
                                  fun = function(y){mean(y, trim = trim)}) +
            ggplot2::stat_summary(ggplot2::aes(group = consistency.metric),
                                  geom = "point",
                                  fun = function(y){mean(y, trim = trim)}) +
            ggplot2::stat_summary(ggplot2::aes(group = consistency.metric),
                                  geom = "errorbar",
                                  fun.data = function(y){trim_mean_se::trim_mean_se(y, trim = trim)},  # Use mean and standard error
                                  width = 0.05) +  # Adjust the width of error bars
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

