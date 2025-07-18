# compute gene markers by:
# 1- Highly variable genes (HVG)
# 2- DEG between clusters (DEG)
# 3- Others...

# This will be added to slot gene.list of scTypeEval object
# Function to compute HVGs for a single matrix
compute_hvg <- function(mat,
                        ngenes = 500,
                        margin = 1L) {
   if (margin == 1L) {
      gene_mean <- Matrix::rowMeans(mat)
      gene_variance <- Matrix::rowSums((mat - gene_mean)^2) / (ncol(mat) - 1)
   } else if (margin == 2L) {
      gene_mean <- Matrix::colMeans(mat)
      gene_variance <- Matrix::colSums((mat - gene_mean)^2) / (nrow(mat) - 1)
   } else {
      stop("Invalid margin. Use 1 for rows or 2 for columns.")
   }
   
   # Coefficient of variation
   cv <- gene_variance / gene_mean
   cv[is.na(cv) | is.infinite(cv)] <- 0  # Handle edge cases
   
   # Rank genes by CV and select the top ngenes
   ranked_genes <- order(cv, decreasing = TRUE)[1:ngenes]
   if (margin == 1L) {
      hvg <- rownames(mat)[ranked_genes]
   } else {
      hvg <- colnames(mat)[ranked_genes]
   }
   return(hvg)
}

get.HVG <- function(norm.mat, 
                    sample = NULL,
                    ngenes = 500, 
                    margin = 1L,
                    bparam = BiocParallel::SerialParam()) {
   # If sample is NULL, compute HVGs for the entire matrix
   if (is.null(sample)) {
      top_hvgs <- compute_hvg(mat = norm.mat,
                              ngenes = ngenes,
                              margin = margin)
   } else {
      
      # If sample is not NULL, compute HVGs for each sample
      hvg_per_sample <- BiocParallel::bplapply(unique(sample),
                                               BPPARAM = bparam,
                                               function(s) {
                                                  sample_mat <- norm.mat[, sample == s, drop = FALSE]
                                                  compute_hvg(mat = sample_mat,
                                                              ngenes = ngenes * 2, # double ngenes to get final list desired
                                                              margin = margin)
                                               })
      
      # Create a table of all HVGs and their frequencies across samples
      all_hvgs <- unique(unlist(hvg_per_sample))
      
      # Calculate frequency of each HVG across samples
      freq_table <- table(unlist(hvg_per_sample))
      hvg_freq <- freq_table[all_hvgs]
      
      # Compute rank within each sample
      rank_matrix <- sapply(hvg_per_sample, function(hvg) {
         ranks <- match(all_hvgs, hvg)  # Rank within this sample
         ranks[is.na(ranks)] <- ngenes * 2     # Set non-present genes to max
         ranks
      })
      
      # Calculate median rank across samples for each gene
      med_rank <- apply(rank_matrix, 1, function(x) median(x, na.rm = TRUE))
      
      # Rank genes first by frequency, then by median rank
      rank_order <- order(-hvg_freq, med_rank)
      top_hvgs <- all_hvgs[rank_order][1:ngenes]

   }
   
   return(top_hvgs)

}


get.GeneVar <- function(norm.mat,
                        sample = NULL,
                        ngenes = 500, 
                        equiweight = TRUE,
                        bparam = BiocParallel::SerialParam(),
                        ...){
   
   var <- scran::modelGeneVar(x = norm.mat,
                              block = sample,
                              BPPARAM = bparam,
                              equiweight = equiweight,
                              ...)
   
   hgv <- scran::getTopHVGs(var,
                            n=ngenes)
      
   return(hgv)
   
}



# get markers using scran findMarkers
get.DEG <- function(mat, # normalized gene expression matrix!,
                    ident,
                    block = NULL,
                    ngenes.celltype = 50,
                    test.type = "t",
                    ncores = 1,
                    bparam = NULL,
                    progressbar = TRUE,
                    min.prop = 0.6,
                    unlist = TRUE, # wheter to return a vector of genes (TRUE) or a list per ident (FALSE)
                    ...){
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   de <- scran::findMarkers(x = mat,
                            groups = ident,
                            block = block,
                            pval.type = "some",
                            BPPARAM = param,
                            min.prop = min.prop,
                            test.type = test.type,
                            ...)
   
   markers <- lapply(de,
                function(d){
                   d <- d |> 
                      as.data.frame() |>
                      dplyr::filter(FDR < 0.05) |>
                      dplyr::arrange(dplyr::desc(summary.logFC)) |>
                      head(ngenes.celltype)
                   
                   return(rownames(d))
                })
   if(unlist){
      markers <- markers |>
         unlist() |>
         unique()
   }

   return(markers)
   
}



