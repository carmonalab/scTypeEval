mutual_method <-  c("Score", "Match")
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
      logfc_vec <- do.call(cbind, lapply(others, function(other) {
         log2((medians[, cl] + pseudo_count) /
                 (medians[, other] + pseudo_count))
      }))
      colnames(logfc_vec) <- others
      
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
                           ref.ident,
                           min.markers = 10,
                           bparam = BiocParallel::SerialParam()){
   
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
                           ref.ident,
                           bparam = BiocParallel::SerialParam()){
   
   # SingleR prediction
   pred <- SingleR::SingleR(
      test = test,
      ref = ref,
      labels = ref.ident,
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

score.tidy <- function(predi,
                       .true){
   pred <- predi |>
      dplyr::select(dplyr::starts_with("score")) |>  
      tibble::rownames_to_column("true") |>
      tidyr::pivot_longer(-c(true),
                          names_to = "celltype",
                          values_to = "score") |>
      dplyr::mutate(celltype = gsub("scores[.]", "", celltype),
                    celltype = names(.true)[match(celltype, .true)]) 
   
   return(pred)
}

.score <- function(pred1, true1,
                   pred2, true2){
   
   pred1 <- score.tidy(pred1, true2)
   pred2 <- score.tidy(pred2, true1)
   names(pred2)[1:2] <- c("celltype", "true")
   
   pred <- dplyr::inner_join(pred1, pred2,
                             by = c("true", "celltype")) |>
      dplyr::mutate(
         score0 = score.x * score.y,
         score = 1 - score0
      ) |>
      dplyr::select(i = true, j = celltype, score)
   
   return(pred)
}

match.tidy <- function(predi, true1, true2){
   # get cell types contained in both samples
   common_celltypes <- intersect(true1, true2)
   
   pred <- predi |>
      dplyr::select("label") |>
      dplyr::mutate(true = true1,
                    # keep as 0.5 when actual match was not possible
                    score = ifelse(true %in% common_celltypes,
                                   0, 0.5)
      ) |>
      tibble::rownames_to_column("id") |>
      dplyr::select(true, label, id, score)
   
   return(pred)
}


.match <- function(pred1, pred2,
                   true1, true2){
   
   pred1 <- match.tidy(pred1, true1 = true1, true2 = true2)
   pred2 <- match.tidy(pred2, true1 = true2, true2 = true1)
   names(pred2)[1:2] <- c("label", "true")
   
   pred <- dplyr::inner_join(pred1, pred2, by = c("true", "label")) |>
      dplyr::rowwise() |>
      dplyr::mutate(score = max(score.x, score.y)) |>
      dplyr::select(i = id.x, j = id.y, score)
   
   return(pred)
}




RecipClassif <- function(mat,
                              ident,
                              sample,
                              group,
                              method = "Score",
                              classifier = "SingleR",
                              ncores = 1,
                              bparam = NULL,
                              progressbar = TRUE,
                              verbose = TRUE){
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   # split by sample
   mat.split <- split_matrix(mat = mat,
                             ident = ident,
                             sample = sample,
                             bparam = param)
   n <- length(mat.split)
   combs <- split(expand.grid(1:n, 1:n), seq_len(n^2))
   
   # set classifier
   classifier <- switch(classifier,
                        "Spearman_correlation" = spear_classify,
                        "SingleR" = singleR.helper,
                        stop(classifier, " classifier not supported")
   )
   
   if(verbose){message("Computing pairwise BesHiT similarity... \n")}
   df.tmp <- BiocParallel::bplapply(combs,
                                    BPPARAM = param,
                                    function(i) {
                                       i <- unlist(i)
                                       a <- i[1]
                                       b <- i[2]
                                       
                                       # Check that both reference and query are not NULL
                                       if (is.null(mat.split[[a]]) || is.null(mat.split[[b]])) {
                                          return(NULL)
                                       }
                                       
                                       # Check if each matrix has at least two unique cell types
                                       nc <- lapply(mat.split[c(a, b)], function(x) length(unique(x@ident)))
                                       if (any(nc < 2)) {
                                          return(NULL)
                                       }
                                       
                                       # Perform classification
                                       pred1 <- classifier(mat.split[[a]]@matrix,
                                                           mat.split[[b]]@matrix,
                                                           ref.ident = mat.split[[b]]@ident)
                                       pred2 <- classifier(mat.split[[b]]@matrix,
                                                           mat.split[[a]]@matrix,
                                                           ref.ident = mat.split[[a]]@ident)
                                       
                                       true1 <- mat.split[[a]]@ident
                                       true2 <- mat.split[[b]]@ident
                                       
                                       # Evaluate metrics
                                       meth <- tolower(method)
                                       results <- switch(meth,
                                                         "score" = .score(pred1 = pred1,
                                                                          pred2 = pred2,
                                                                          true1 = true1,
                                                                          true2 = true2),
                                                         "match" = .match(pred1 = pred1,
                                                                          pred2 = pred2,
                                                                          true1 = true1,
                                                                          true2 = true2),
                                                         stop(meth, " is not a supported Reciprocal Classficiation method."))
                                       
                                       lst <- apply(results, 1, as.list)
                                       
                                       # return a data frame with one column with cell type and one column with score
                                       return(lst)
                                    })
   
   # Filter out NULLs
   df.tmp <- Filter(Negate(is.null), df.tmp)
   df.tmp <- unlist(df.tmp, recursive = F)
   
   # ---- Build pairwise matrix with lists
   n <- length(group)
   # Initialize distance matrix
   dist_mat <- matrix(1, n, n)
   rownames(dist_mat) <- colnames(dist_mat) <- group
   
   
   # Fill in symmetric matrix
   for (a in seq_along(df.tmp)) {
      res <- df.tmp[[a]]
      i <- res$i
      j <- res$j
      d <- res$score
      dist_mat[i, j] <- dist_mat[j, i] <- as.numeric(d)
   }
   
   
   return(as.dist(dist_mat))
   
}
