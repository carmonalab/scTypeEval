mutual_method <-  c("Mutual.Score", "Mutual.Match")
classifiers <- c("Spearman_correlation", "SingleR")

get.DEG_logfc <- function(mat,
                          ident,
                          ngenes.celltype = NULL,
                          min.logfc = 1,
                          pseudo_count = 1e-6) {
   
   clusters <- levels(ident)
   
   if(is.null(ngenes.celltype)){
      ngenes.celltype <- floor(500 * (2/3)^log2(length(clusters)))
   }
   
   # Normalize 
   
   # Initialize a list to store medians per cluster
   medians_list <- list()
   # Compute row medians for each cluster
   for (cl in clusters) {
      cols <- which(ident == cl)
      dense_mat <- as.matrix(mat[, cols, drop = FALSE])
      medians_list[[cl]] <- apply(dense_mat, 1, median)
   }
   
   # Combine medians into a matrix
   medians <- do.call(cbind, medians_list)
   
   # Initialize marker list
   markers <- vector("list", length(clusters))
   names(markers) <- clusters
   
   # For each cluster, find top genes based on pairwise log2FC
   for (cl in clusters) {
      others <- setdiff(clusters, cl)
      
      # Store the maximum logFC against all other clusters
      logfc_vec <- sapply(others, function(other) {
         log2((medians[, cl] + pseudo_count) / 
                 (medians[, other] + pseudo_count))
      })
      
      top_de <- lapply(others, function(x){
         m <- logfc_vec[,x]
         # keep only FC > min.logfc
         m <- m[m > min.logfc]
         # filter lowly expressed genes
         # cols <- which(ident == x)
         # emat <- mat[, cols, drop = FALSE]
         # kg <- rownames(emat[Matrix::rowSums(emat) >= length(cols), ,drop = F])
         # m <- m[names(m) %in% kg]
         # sort by FC
         m <- sort(m, decreasing = TRUE)
         nge <- min(ngenes.celltype, length(m))
         names(m[1:nge])
      }) |>
         unlist() |>
         unique()
      
      markers[[cl]] <- top_de
   }
   return(markers)
}


# function to classify cell type based on query-cells vs ref-cells spearman correlation
# test and ref must be sparse matrices
spear_classify <- function(test_mat,
                           ref_mat,
                           min.markers = 10,
                           bparam = BiocParallel::SerialParam()){
   test <- test_mat@matrix |>
      Normalize_data()
   ref <- ref_mat@matrix |>
      Normalize_data()
   ref.ident <- ref_mat@ident
   
   # Step 0: ensure overlapping genes
   common_genes <- intersect(rownames(test), rownames(ref))
   test <- test[common_genes, , drop = FALSE]
   ref <- ref[common_genes, , drop = FALSE]
   
   # Step 1: get DEG markers based only on median log2FC
   markers <- get.DEG_logfc(ref,
                            ident = ref.ident)
   
   # Step 2: group reference by identity
   ref_groups <- split(seq_len(ncol(ref)), ref.ident)
   
   # Step 3: compute scores per identity
   scores_list <- 
      BiocParallel::bplapply(names(ref_groups),
                             BPPARAM = bparam,
                             function(ident_name) {
                                marker_genes <- intersect(markers[[ident_name]], rownames(test))
                                
                                if (length(marker_genes) < min.markers) {
                                   return(rep(NA, ncol(test)))  # skip if too few markers
                                }
                                
                                ref_mat <- ref[marker_genes, ref_groups[[ident_name]], drop = FALSE]
                                test_mat <- test[marker_genes, , drop = FALSE]
                                
                                apply(test_mat, 2, function(test_col) {
                                   quantile(apply(ref_mat, 2, function(ref_col) {
                                      # get percentil 80th of correlations
                                      cor(ref_col, test_col, method = "spearman")
                                   }), probs = 0.8, na.rm = TRUE)
                                })
                             })
   
   # Step 4: build score dataframe
   score_df <- as.data.frame(do.call(cbind, scores_list))
   colnames(score_df) <- paste0("scores.", names(ref_groups))
   rownames(score_df) <- colnames(test)
   
   # Step 5: predicted label (handling NA rows)
   score_df$label <- apply(score_df, 1, function(x) {
      if (all(is.na(x))) {
         return(NA)
      } else {
         names(ref_groups)[which.max(x)]
      }
   })
   
   return(score_df)
}



# code to perform a best-hit classifier based on singleR
singleR.helper <- function(test,
                           ref,
                           bparam = BiocParallel::SerialParam()){
   # transform to SCE
   test <- get.SCE(test)
   ref <- get.SCE(ref)
   # SingleR prediction
   pred <- SingleR::SingleR(
      test = test,
      ref = ref,
      labels = ref$ident,
      de.method = "classic",
      BPPARAM = bparam
   ) |> 
      as.data.frame() |>
      # condition for when singleR return no results (no label)
      dplyr::mutate(label = ifelse(is.na(pruned.labels),
                                   "_nores_",
                                   pruned.labels
      )
      )
   return(pred)
}

score.tidy <- function(pred,
                       .true,
                       filter = TRUE){
   pred <- pred |>
      dplyr::select(dplyr::starts_with("score")) |>  
      dplyr::mutate(true = .true) |>
      tidyr::pivot_longer(-c(true),
                          names_to = "celltype",
                          values_to = "score") 
   if (filter) {
      pred <- pred |> 
         dplyr::filter(true == gsub("scores[.]", "", celltype)) |> 
         dplyr::group_by(true) |> 
         dplyr::summarize(score = mean(score))
   }
   
   return(pred)
}

.score <- function(pred1, true1,
                   pred2, true2){
   
   pred1 <- score.tidy(pred1, true1)
   pred2 <- score.tidy(pred2, true2)
   
   pred <- dplyr::inner_join(pred1, pred2, by = "true") |>
      dplyr::mutate(score = score.x * score.y) |>
      dplyr::select(true, score) |>
      dplyr::rename("celltype" = "true")
   
   return(pred)
}

match.tidy <- function(pred, .true){
   pred <- pred |>
      dplyr::select("label") |>
      dplyr::mutate(true = .true,
                    score = ifelse(true == label,
                                   1, 0
                    )
      )
   
   return(pred)
}


.match <- function(pred1, pred2,
                   true1, true2){
   
   pred1 <- match.tidy(pred1, true1)
   pred2 <- match.tidy(pred2, true2)
   
   pred <- dplyr::inner_join(pred1, pred2, by = "true") |>
      dplyr::group_by(true) |>
      dplyr::mutate(score = score.x * score.y) |>
      dplyr::select(true, score) |>
      dplyr::rename("celltype" = "true")
   
   return(pred)
}

get.MutualMatrix <- function(mat,
                             ident,
                             sample,
                             data.type = "sc",
                             min.cells = 10,
                             min.samples = 5,
                             sep = "_",
                             bparam = BiocParallel::SerialParam()){
   # Combined grouping
   groups <- factor(paste(sample, ident, sep = sep))
   
   valid_indices <- valid.indices(groups = groups,
                                  ident = ident,
                                  sample = sample,
                                  min.samples = min.samples,
                                  min.cells = min.cells,
                                  sep = sep)
   
   mat <- mat[, valid_indices, drop = FALSE]
   
   mats <- split_matrix(mat,
                        ident = ident[valid_indices],
                        sample = sample[valid_indices],
                        min.cells = min.cells,
                        bparam = bparam)
   
   mat.split <- lapply(mats,
                       function(mat){
                          get.matrix(mat@matrix,
                                     data.type = data.type,
                                     ident = mat@ident,
                                     sample = NULL, # now is pseudobulk per cell type per invidividual sample
                                     min.cells = min.cells,
                                     bparam = bparam)
                       })
   
   return(mat.split)
}

get.SCE <-  function(m){
   sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = m@matrix),
      colData = data.frame(ident = m@ident)
   )
   #sce <- sce[,Matrix::colSums(SingleCellExperiment::counts(sce)) > 0] # Remove libraries with no counts.
   # this is taken account in the processing steps
   sce <- scuttle::logNormCounts(sce)
   return(sce)
}

bestHit <- function(mat,
                    ident,
                    ident_GroundTruth = NULL,
                    sample,
                    method = "Mutual.Score",
                    classifier = "Spearman_correlation",
                    data.type = "sc",
                    min.cells = 10,
                    min.samples = 5,
                    sep = "_",
                    ncores = 1,
                    bparam = NULL,
                    progressbar = TRUE){
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   mat.split <- get.MutualMatrix(mat = mat,
                                 ident = ident,
                                 sample = sample,
                                 data.type = data.type,
                                 min.cells = min.cells,
                                 min.samples = min.samples,
                                 sep = sep,
                                 bparam = param)
   
   if(is.null(ident_GroundTruth)){
      mat.splitGT <- mat.split
   } else {
      mat.splitGT <- get.MutualMatrix(mat = mat,
                                      ident = ident_GroundTruth,
                                      sample = sample,
                                      data.type = data.type,
                                      min.cells = min.cells,
                                      min.samples = min.samples,
                                      sep = sep,
                                      bparam = param)
   }
   
   
   names_list <- names(mat.split)
   idx <- expand.grid(i = seq_along(names_list),
                      j = seq_along(names_list)) |>
      filter(i < j)
   
   combis <- data.frame(Var1 = names_list[idx$i],
                        Var2 = names_list[idx$j])
   
   # set classifier
   classifier <- switch(classifier,
                        "Spearman_correlation" = spear_classify,
                        "SingleR" = singleR.helper,
                        stop(classifier, " classifier not supported")
   )
   
   df.tmp <- list()
   
   df.tmp <- BiocParallel::bplapply(seq_len(nrow(combis)),
                                    BPPARAM = param,
                                    function(i) {
      a <- combis[i, 1]
      b <- combis[i, 2]
      
      # Check that both reference and query are not NULL
      if (is.null(mat.split[[a]]) || is.null(mat.split[[b]]) ||
          is.null(mat.splitGT[[a]]) || is.null(mat.splitGT[[b]])) {
         return(NULL)
      }
      
      # Check if each matrix has at least two unique cell types
      nc <- lapply(mat.splitGT[c(a, b)], function(x) length(unique(x@ident)))
      if (any(nc < 2)) {
         return(NULL)
      }
      
      # Perform classification
      pred1 <- classifier(mat.split[[a]], mat.splitGT[[b]])
      pred2 <- classifier(mat.split[[b]], mat.splitGT[[a]])
      
      true1 <- mat.split[[a]]@ident
      true2 <- mat.split[[b]]@ident
      
      # Evaluate metrics
      results <- lapply(method, function(meth) {
         switch(meth,
                "Mutual.Score" = .score(pred1 = pred1,
                                        pred2 = pred2,
                                        true1 = true1,
                                        true2 = true2),
                "Mutual.Match" = .match(pred1 = pred1,
                                        pred2 = pred2,
                                        true1 = true1,
                                        true2 = true2),
                stop(meth, " is not a supported Mutual Besthit method."))
      })
      names(results) <- method
      return(results)
   })
   
   # Filter out NULLs
   df.tmp <- Filter(Negate(is.null), df.tmp)
   
   # join results by method
   resli <- lapply(method, function(m)
   {
      res <- do.call(rbind, lapply(df.tmp, function(x) x[[m]]))
      # summarize results per cell type
      res <- res |>
         dplyr::group_by(celltype) |>
         dplyr::summarize(score = mean(score, na.rm = T))
      
      if(m == "Mutual.Match"){
         # normalize for BestHiT match
         nident <- sapply(mat.split, function(mat) {length(mat@ident)})
         # Generate all pairwise combinations of the samples
         sample_combinations <- combn(nident, 2)
         # Calculate pairwise probabilities (1 / (cell types in sample i * cell types in sample j))
         probabilities <- apply(sample_combinations, 2, function(pair) {
            1 / (pair[1] * pair[2])})
         # Average all pairwise probabilities together
         min.val <- mean(probabilities) 
         res <- res |>
            dplyr::mutate(score = minmax_norm(score,
                                              min_value = min.val,
                                              max_value = 1))
      }
      # extract vector with output
      res <-  res |>
         dplyr::pull(score, name = celltype)
      
      return(res)
   }
   )
   names(resli) <- method
   
   return(resli)
   
}
