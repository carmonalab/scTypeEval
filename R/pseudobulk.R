# script to pseudobulking cells:

# get either regular pseudobulk or 1 vs all pseudobulk
get_pseudobulk <- function(mat, 
                           ident = NULL, 
                           sample = NULL, 
                           min.samples = 5,  # Minimum number of samples with > 5 cells for each ident
                           min.cells = 10,   # Minimum cells per pseudobulk to be included
                           ncores = 1, 
                           sep = "_") {  
   
   if (is.null(sample)) {
      # Single grouping
      groups <- factor(ident)
   } else {
      # Combined grouping
      groups <- factor(paste(sample, ident, sep = sep))
   }
   
   group_levels <- levels(groups)
   group_indices <- as.integer(groups)
   
   # Count the number of cells per group
   cell_counts <- table(groups)
   
   # Filter groups with at least `min.cells` cells
   valid_groups <- names(cell_counts[cell_counts >= min.cells])
   
   # Check that each ident has at least `min.samples` samples with more than 5 cells
   if (!is.null(sample)) {
      # Split by ident and check sample distribution
      ident_samples_count <- sapply(unique(ident), function(cell_type) {
         sample_counts <- table(sample[ident == cell_type])
         sum(sample_counts >= 5)  # Count samples with at least 5 cells for this cell type
      })
      
      valid_idents <- names(ident_samples_count[ident_samples_count >= min.samples])
      
      # Filter valid groups based on the number of samples with at least 5 cells
      valid_groups <- valid_groups[sapply(valid_groups, function(g) {
         ident_name <- strsplit(g, sep)[[1]][2]
         ident_name %in% valid_idents
      })]
   }
   
   # Filter the grouping indices to keep only valid groups
   valid_indices <- which(groups %in% valid_groups)
   if (length(valid_indices) == 0) {
      stop("No valid pseudobulk groups meet the specified thresholds.")
   }
   
   groups <- groups[valid_indices]
   mat <- mat[, valid_indices, drop = FALSE]
   
   group_levels <- levels(groups)
   group_indices <- as.integer(groups)
   
   # Create a sparse matrix for grouping, ensuring groups are columns
   grouping_matrix <- Matrix::sparseMatrix(
      i = seq_along(group_indices),       # Set rows as indices for each sample
      j = group_indices,                  # Set columns as the group indices
      x = 1,                              # Weights are 1 (indicating inclusion in the group)
      dims = c(ncol(mat), length(group_levels))  # Dimensions: genes (rows) x groups (columns)
   )
   
   # Multiply the grouping matrix with the count matrix to get the summed matrix
   summed_matrix <- grouping_matrix %*% mat
   
   # Set row names for clarity (genes as rows, groups as columns)
   colnames(summed_matrix) <- group_levels
   
   return(summed_matrix)
}


# sketching 1 vs all at single-cel level

get_PCA <- function(norm.mat,
                      hgv,
                      ndim = 30){
   # compute PCA
   pca <- prcomp(t(norm.mat[hgv,]),
                 rank. = ndim)
   return(pca)
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
