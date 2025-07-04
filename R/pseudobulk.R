# script to pseudobulking cells:

data_type <- c("sc", "pseudobulk",
               "pseudobulk_1vsall",
               "GloScope")


valid.indices <- function(groups,
                          ident,
                          sample,
                          min.cells = 10,
                          min.samples = 5,
                          sep = "_"){
   
   
   # Count the number of cells per group
   cell_counts <- table(groups)
   
   # Filter groups with at least `min.cells` cells
   valid_groups <- names(cell_counts[cell_counts >= min.cells])
   
   #update ident
   nident <- ident[groups %in% valid_groups]
   
   # Check that each ident has at least `min.samples` samples with more than `min.cells` cells
   if (!is.null(sample)) {
      # update samples
      nsample <- sample[groups %in% valid_groups]
      # Split by ident and check sample distribution
      ident_samples_count <- sapply(unique(nident), function(cell_type) {
         sample_counts <- table(nsample[nident == cell_type])
         sum(sample_counts >= min.samples)  # Count samples with at least `min.samples` cells for this cell type
      })
      
      names(ident_samples_count) <- unique(nident)
      
      valid_idents <- names(ident_samples_count[ident_samples_count >= min.samples])
      
      # Filter valid groups based on the number of samples with at least `min.cells` cells
      valid_groups <- valid_groups[sapply(valid_groups, function(g) {
         ident_name <- strsplit(g, sep)[[1]][2]
         ident_name %in% valid_idents
      })]
   }
   
   # Filter the grouping indices to keep only valid groups
   valid_indices <- which(groups %in% valid_groups)
   if (length(valid_indices) == 0) {
      stop("No valid sample-celltype groups meet the specified thresholds.")
   }
   
   return(valid_indices)
}


# get either regular pseudobulk or 1 vs all pseudobulk
get_pseudobulk <- function(mat, 
                           ident, 
                           sample = NULL, 
                           min.samples = 5,  # Minimum number of samples with > 5 cells for each ident
                           min.cells = 10,   # Minimum cells per pseudobulk to be included 
                           sep = "_") {  
   
   if (is.null(sample)) {
      # Single grouping
      groups <- factor(ident)
   } else {
      # Combined grouping
      groups <- factor(paste(sample, ident, sep = sep))
   }
   
   valid_indices <- valid.indices(groups = groups,
                                  ident = ident,
                                  sample = sample,
                                  min.samples = min.samples,
                                  min.cells = min.cells,
                                  sep = sep)

   groups <- factor(groups[valid_indices])
   mat <- mat[, valid_indices, drop = FALSE]
   
   group_levels <- levels(groups)
   group_indices <- as.integer(groups)
   
   grouping_matrix <- Matrix::sparseMatrix(
      i = seq_along(group_indices),           # Cell indices (columns)
      j = group_indices,                      # Group indices (rows)
      x = 1,                                  # Weights are 1 (indicating inclusion in the group)
      dims = c(ncol(mat), length(group_levels))  # Cells (columns) x Groups (rows)
   )
   
   # Multiply the grouping matrix with the count matrix to get the summed matrix
   summed_matrix <- mat %*% grouping_matrix
   
   # Set row names for clarity (genes as rows, groups as columns)
   colnames(summed_matrix) <- group_levels
   
   # Filter out pseudobulks (columns) with total count of 0
   keep_cols <- Matrix::colSums(summed_matrix) > 0
   summed_matrix <- summed_matrix[, keep_cols, drop = FALSE]
   group_levels <- group_levels[keep_cols]
   
   # retrieve new idents
   if(is.null(sample)){
      new.idents <- group_levels
   } else {
      new.idents <- sapply(group_levels, function(x){strsplit(x, sep)[[1]][2]})
   }
   
   ret <- new("Mat_ident",
              matrix = summed_matrix,
              ident = factor(new.idents))
   
   return(ret)
}


get.pseudobulk1vsAll <- function(mat, 
                                 ident = NULL, 
                                 sample = NULL, 
                                 min.samples = 5,  # Minimum number of samples with > 5 cells for each ident
                                 min.cells = 10, 
                                 ncores = NULL,
                                 bparam = NULL,
                                 sep = "_"){
   # check if pseudobulk 1 vs all is possible
   if(is.null(sample)){
      stop("1vsAll pseudobulk only possible for multiple samples.
           Please provide `sample` parameter.")
   }
   
   groups <- factor(paste(sample, ident, sep = sep))
   names(groups) <- sample
   
   
   valid_indices <- valid.indices(groups = groups,
                                  ident = ident,
                                  sample = sample,
                                  min.samples = min.samples,
                                  min.cells = min.cells,
                                  sep = sep)
   ngroups <- groups[valid_indices]
   nident <- sapply(as.character(ngroups),
                    function(x)strsplit(x, sep)[[1]][2]) |> unique()
   
   nsamples <- lapply(nident,
                      function(x){
                         g <- sapply(as.character(ngroups),
                                     function(x){strsplit(x, sep)[[1]][2]})
                         names(g) <- sapply(as.character(ngroups),
                                     function(x){strsplit(x, sep)[[1]][1]})
                         g <- g[g == x]
                         r <- unique(names(g))
                         return(r)
                      })
   names(nsamples) <- nident
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   mat.list <- BiocParallel::bplapply(nident,
                                      BPPARAM = param,
                                      function(i){
                                         vsample <- nsamples[[i]]
                                         ind <- names(groups) %in% vsample
                                         
                                         nmat <- mat[,ind]
                                         ident <- ident[ind]
                                         sample <- sample[ind]
                                         
                                         ident <- as.character(ident)
                                         ident[ident != i] <- "psblk"
                                         ident <- as.factor(ident)
                                         
                                         mat <- get_pseudobulk(nmat, 
                                                               ident = ident, 
                                                               sample = sample,
                                                               min.samples = min.samples,
                                                               min.cells = min.cells)
                                         
                                         return(mat)
                                      })
   
   names(mat.list) <- nident
   
   return(mat.list)
   
}

valid.sc <- function(mat,
                     ident,
                     min.cells = 10){
   # Count the number of cells per group
   cell_counts <- table(ident)
   # Filter groups with at least `min.cells` cells
   valid_groups <- names(cell_counts[cell_counts >= min.cells])
   keep_group <- ident %in% valid_groups
   
   # Remove cells with total count of 0
   keep_nonzero <- Matrix::colSums(mat) > 0
   
   # Combine both filters
   keep <- keep_group & keep_nonzero
   
   # Subset matrix and ident vector
   mat <- mat[, keep, drop = FALSE]
   new.ident <- ident[keep]
   
   ret <- new("Mat_ident",
              matrix = mat,
              ident = factor(new.ident))
   
   return(ret)
}


get.sc <- function(mat,
                   ident,
                   sample = NULL,
                   min.cells = 10,
                   bparam = BiocParallel::SerialParam()){
   
   if(is.null(sample)){
      ret <- valid.sc(mat = mat,
                      ident = ident,
                      min.cells = min.cells)
   } else {
      ret <- split_matrix(mat = mat,
                          ident = ident,
                          sample = sample,
                          bparam = bparam)
   }
   
   return(ret) 
   
}


# Keep the filtering of cell type with low proportion or in few samples
get.GloScope <- function(mat, 
                         ident, 
                         sample = NULL, 
                         min.samples = 5,  # Minimum number of samples with > 5 cells for each ident
                         min.cells = 10,   # Minimum cells per pseudobulk to be included 
                         sep = "_") {  
   
   if (is.null(sample)) {
      # Single grouping
      groups <- factor(ident)
   } else {
      # Combined grouping
      groups <- factor(paste(sample, ident, sep = sep))
   }
   
   valid_indices <- valid.indices(groups = groups,
                                  ident = ident,
                                  sample = sample,
                                  min.samples = min.samples,
                                  min.cells = min.cells,
                                  sep = sep)
   
   # Remove cells with total count of 0
   keep_nonzero <- Matrix::colSums(mat) > 0
   
   vl <- colnames(mat) %in% colnames(mat)[valid_indices]
   
   keep <- vl & keep_nonzero
   
   groups <- factor(groups[keep])
   mat <- mat[, keep, drop = FALSE]
   
   
   ret <- new("Mat_ident",
              matrix = mat,
              ident = factor(groups))
   
   return(ret)
}


# get matrix according to the data type
get.matrix <- function(matrix,
                       data.type,
                       ident,
                       sample,
                       min.samples = 5,  # Minimum number of samples with > 5 cells for each ident
                       min.cells = 10, # Minimum cells per pseudobulk to be included 
                       bparam = BiocParallel::SerialParam()){
   
   
   mat <- switch(data.type,
                 "sc" = get.sc(mat = matrix,
                               ident = ident,
                               sample = sample,
                               min.cells = min.cells,
                               bparam = bparam),
                 "pseudobulk" = get_pseudobulk(mat = matrix, 
                                               ident = ident, 
                                               sample = sample,
                                               min.samples = min.samples,
                                               min.cells = min.cells),
                 # return a list of matrix, one for each ident
                 "pseudobulk_1vsall" = get.pseudobulk1vsAll(mat = matrix, 
                                                            ident = ident, 
                                                            sample = sample,
                                                            min.samples = min.samples,
                                                            min.cells = min.cells,
                                                            bparam = bparam),
                 "GloScope" = get.GloScope(mat = matrix,
                                           ident = ident,
                                           sample = sample,
                                           min.cells = min.cells,
                                           min.samples = min.samples),
                 stop(data.type, " is not a supported data.type.")
                 )
                 
   return(mat)
}



# get_intra_SvB <- function(object,
#                           colLabel,
#                           celltype,
#                           min.cells = 10,
#                           nfeatures = 1000,
#                           seed = 22,
#                           ncores = 1) {
#    
#    # set paralelization
#    param <- set_parallel_params(ncores = ncores, progressbar = F)
#    
#    ct.mat <- object[, object@meta.data[[colLabel]] == celltype][[assay]]$counts
#    colnames(ct.mat) <- gsub("_", "-", colnames(ct.mat)) %>%
#       paste(., celltype, sep = "_")
#    
#    ncells <- ncol(ct.mat)
#    
#    # adjust if their proportions is very discrepant
#    seu.psblk <- object[, object@meta.data[[colLabel]] != celltype]
#    
#    # find variable featurs
#    seu.psblk <- FindVariableFeatures(seu.psblk,
#                                      selection.method = "vst",
#                                      nfeatures = nfeatures,
#                                      verbose = F)
#    
#    # produce sketching and merge recursively
#    ps.mat <- bplapply(1:ncells,
#                       BPPARAM = param,
#                       function(n){
#                          # sketching
#                          
#                          myseed <- seed + n
#                          suppressMessages(
#                             {
#                                ds <- Seurat::SketchData(seu.psblk,
#                                                         assay = assay,
#                                                         ncells = min.cells,
#                                                         seed = myseed,
#                                                         verbose = F)
#                                ds.mat <- ds[["sketch"]]$counts
#                                ds.mat <- Matrix::rowSums(ds.mat)
#                             })
#                          
#                          return(ds.mat)
#                       })
#    
#    # join Sketch pseudobulk
#    ps.mat <- do.call(cbind, ps.mat)
#    
#    colnames(ps.mat) <- paste0(paste0("r", 1:ncells), "_psblk")
#    
#    # join with actual celltypes
#    fi <- cbind(ct.mat, ps.mat)
#    
#    return(fi)
#    
#    
# }
