# script to pseudobulking cells:

aggregation_types <- c("single-cell", "pseudobulk")


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

   
   # retrieve new idents, group and sample
   new.group <- factor(group_levels)
   new.ident <- sapply(group_levels, function(x){strsplit(x, sep)[[1]][2]})
   new.ident <- factor(new.ident)
   
   if(!is.null(sample)){
      new.sample <- sapply(group_levels, function(x){strsplit(x, sep)[[1]][1]})
      new.sample <- factor(new.sample)
   } else {
      new.sample <- NULL
   }
   
   ret <- new("Mat_ident",
              matrix = summed_matrix,
              group = new.group,
              ident = new.ident,
              sample = new.sample)
   
   return(ret)
}



# Keep the filtering of cell type with low proportion or in few samples
get_singlecell <- function(mat, 
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
   
   
   mat <- mat[, valid_indices, drop = FALSE]
   
   new.group <- factor(groups[valid_indices])
   new.ident <- factor(ident[valid_indices])
   if(!is.null(sample)){
      new.sample <- factor(sample[valid_indices])
   } else {
      new.sample <- NULL
   }
   
   
   ret <- new("Mat_ident",
              matrix = mat,
              group = new.group,
              ident = new.ident,
              sample = new.sample)
   
   return(ret)
}


# get matrix according to the data type
get.matrix <- function(matrix,
                       aggregation,
                       ident,
                       sample,
                       min.samples = 5,  # Minimum number of samples with > 5 cells for each ident
                       min.cells = 10 # Minimum cells per sample-ident to be included 
                       ){
   
   
   mat <- switch(aggregation,
                 "single-cell" = get_singlecell(mat = matrix,
                                        ident = ident,
                                        sample = sample,
                                        min.samples = min.samples,
                                        min.cells = min.cells
                                        ),
                 "pseudobulk" = get_pseudobulk(mat = matrix, 
                                               ident = ident, 
                                               sample = sample,
                                               min.samples = min.samples,
                                               min.cells = min.cells
                                               ),
                 stop(data.type, " is not a supported aggregation method.")
   )
   
   return(mat)
}


# function to filter emtpy cols or rows in matrix
filter_empty <- function(mat_ident){
   # Extract slots
   mat <- mat_ident@matrix
   group <- mat_ident@group
   ident.name <- names(mat_ident@ident)
   ident <- unlist(mat_ident@ident)
   sample <- mat_ident@sample
   
   # Identify non-zero rows and columns
   keep_rows <- Matrix::rowSums(mat) > 0
   keep_cols <- Matrix::colSums(mat) > 0
   
   # Filter matrix
   mat_filtered <- mat[keep_rows, keep_cols, drop = FALSE]
   
   # Filter metadata
   group_filtered <- factor(group[keep_cols])
   ident_filtered <- factor(ident[keep_cols])
   names(ident_filtered) <- group_filtered
   sample_filtered <- factor(sample[keep_cols])
   names(sample_filtered) <- group_filtered
   
   # Return new Mat_ident object
   ret <- new("Mat_ident",
              matrix = mat_filtered,
              group = group_filtered,
              ident = ident_filtered,
              sample = sample_filtered)
   
   return(ret)
   
}
