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
get.DEG <- function(mat,
                    ident,
                    block = NULL,
                    ngenes.celltype = 50,
                    test.type = "t",
                    ncores = 1,
                    bparam = NULL,
                    progressbar = TRUE,
                    min.prop = 0.6,
                    black.list = NULL,
                    ...){
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   norm.mat <- Normalize_data(mat = mat,
                              method = "Log1p")
   

   norm.mat <- norm.mat[!rownames(norm.mat) %in% black.list,]
   
   de <- scran::findMarkers(x = norm.mat,
                            groups = ident,
                            block = block,
                            pval.type = "some",
                            BPPARAM = param,
                            min.prop = min.prop,
                            ...)
   
   markers <- lapply(de,
                function(d){
                   d <- d |> 
                      as.data.frame() |>
                      dplyr::filter(FDR < 0.05) |>
                      dplyr::arrange(dplyr::desc(summary.logFC)) |>
                      head(ngenes.celltype)
                   
                   return(rownames(d))
                }) |>
      unlist() |>
      unique()

   return(markers)
   
   
}


# gpFISH
get.gpsFISHMarkers <- function(sc_count,
                               ident,
                               rm_unannot = T,
                               ngenes.celltype = 50,
                               total.ngenes = 500,
                               black.list = NULL,
                               seed = 22,
                               ncores = 1){
   # Check if gpsFISH is installed
   if (!requireNamespace("gpsFISH", quietly = TRUE)) {
      stop("The 'gpsFISH' package is required but not installed. Install it with: 
         remotes::install_github('kharchenkolab/gpsFISH') OR
         devtools::install_github('kharchenkolab/gpsFISH')  ")
   }
   
   #https://htmlpreview.github.io/?https://github.com/kharchenkolab/gpsFISH/blob/main/doc/gene_panel_selection.html
   simulation_params <- gpsFISH::simulation_params
   
   # remove black listed genes
   sc_count <- sc_count[!rownames(sc_count) %in% black.list,]
   
   # relative expression
   unique_cluster_label <- as.character(unique(ident))#unique cell type labels
   
   # adapt ident for function
   ccc <- data.frame(cell_name = colnames(sc_count),
                     class_label = ident,
                     row.names = colnames(sc_count))
   
   #cluster-wise relative expression
   relative_prop <-  list()
   relative_prop[["cluster.average"]] <- sapply(unique_cluster_label,
                                                gene_list=rownames(sc_count),
                                                gpsFISH::relative_freq,
                                                count_matrix=sc_count,             
                                                cell_cluster_conversion=ccc)
   #individual cell level relative expression
   relative_prop[["cell.level"]] = t(t(sc_count)/colSums(sc_count))
   
   # Filter genes
   ave_count_ik_all_gene = sapply(unique_cluster_label,
                                  gene_list = rownames(sc_count),
                                  gpsFISH::average_count,
                                  count_matrix = sc_count,               
                                  cell_cluster_conversion = ccc)
   maxexpr <- apply(ave_count_ik_all_gene, 1, max) 
   
   #keep genes with maxexpr >= 1
   gene2keep <- rownames(sc_count)[which(maxexpr>=1)]
   
   if(rm_unannot){
      #remove genes with . in their gene symbol
      gene2remove1 <- grep("\\.", gene2keep)
      #remove lncRNAs (genes with -AS in their gene symbol)
      gene2remove2 <- grep("-AS", gene2keep)
      if (length(union(gene2remove1, gene2remove2))>0){
         gene2keep <- gene2keep[-union(gene2remove1, gene2remove2)]
      }
   }
   
   # Subset the datasets based on genes passing all the filters
   sc_count <- sc_count[gene2keep, ]
   relative_prop$cluster.average <- relative_prop$cluster.average[gene2keep, ]
   relative_prop$cell.level <- relative_prop$cell.level[gene2keep, ]
   
   # Construct weighted penalty matrix based on cell type hierarchy
   cluster_distance <- gpsFISH::cluster_dis(count_table = sc_count,                       
                                            cell_cluster_conversion = ccc,
                                            cluster_metric = "complete",
                                            dist_metric = "correlation",
                                            log.transform = F,                                 
                                            top.var.gene = 1000)
   
   #Then we can construct a raw weighted penalty matrix based on the cell type distance.
   raw_weight_penalty <- as.matrix(cluster_distance$distance_matrix)
   raw_weight_penalty <- raw_weight_penalty/max(raw_weight_penalty) #normalize the value to make sure it is between 0 and 1
   
   weight_penalty <-  gpsFISH::hierarchical_penalty(weight.matrix = raw_weight_penalty,
                                                    cell.type.hierarchy = cell_type_hierarchy,
                                                    reference.resolution = "class",
                                                    current.resolution = "subclass",
                                                    penalty = 2)
   
   #change the diagonal value to 1. This is to make sure that correct predictions will stay unchanged
   diag(weight_penalty) <- 1
   
   # generate gene weight
   #Calculate DEGs for each cell type
   adjust_variance <- preprocess_normalize(sc_count,
                                           n.core= ncores)

   diff_expr=suppressMessages(gpsFISH::diff_gene_cluster(pagoda_object = adjust_variance$pagoda.object,
                                                         cell_cluster_conversion = sc_cluster,
                                                         n.core = ncores))     
   diff_expr_result = diff_expr$diff_result
   
   #statistics for the population
   pop_size = ngenes.celltype #population size: number of gene panels in a population
   panel_size = total.ngenes        #panel size: number of genes in a gene panel
   num.random.panel = 95   #number of panels that we initialize randomly
   num.DE.panel = pop_size-num.random.panel    #number of panels initialized using DEGs
   
   #initialize panels with DEGs
   set.seed(seed)
   initpop.DE <- gpsFISH::initialize_population(pop.size = num.DE.panel,
                                                panel.size = panel_size,
                                                diff_expr_result = diff_expr_result,
                                                diff_metric = "AUC",
                                                diff_metric_cutoff = 0.7,
                                                gene.list = rownames(sc_count),
                                                gene2include = gene2include.symbol)
   #initialize panels with randomly selected genes
   set.seed(seed)
   initpop.random <- gpsFISH::initialize_population_random(pop.size = num.random.panel,
                                                           panel.size = panel_size,
                                                           gene.list = rownames(sc_count),
                                                           gene2include.id = gene2include.id)
   initpop = rbind(initpop.DE, initpop.random)
   
   
   # Gene panel selection
   simulation_params$data2fit = NA
   simulation_params$model.summary = NA
   
   GA = gpsFISH::gpsFISH_optimize(earlyterm = 10,
                                  converge.cutoff = 0.01,
                                  n = dim(sc_count)[1],
                                  k = panel_size,
                                  ngen = 10,
                                  popsize = pop_size,
                                  verbose = 1,
                                  cluster = 1,
                                  initpop = initpop,
                                  method = "NaiveBayes",
                                  metric = "Accuracy",
                                  nCV = 5,
                                  rate = 1,
                                  cluster_size_max = 50,
                                  cluster_size_min = 30,
                                  two_step_sampling_type = c("Subsampling_by_cluster", "Simulation"),
                                  simulation_model = "ZINB",
                                  sample_new_levels = "old_levels",
                                  use_average_cluster_profiles = FALSE,
                                  save.intermediate = FALSE,
                                  full_count_table = as.data.frame(t(sc_count)),
                                  cell_cluster_conversion = sc_cluster,       
                                  relative_prop = relative_prop,
                                  simulation_parameter = simulation_params,
                                  gene2include.id = gene2include.id,
                                  gene.weight = gene.weight,
                                  weight_penalty = weight_penalty
   )
   
   marker_panel = rownames(sc_count)[GA$bestsol]
   
   return(marker_panel)
}
