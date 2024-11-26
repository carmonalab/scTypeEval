# compute gene markers by:
# 1- Highly variable genes (HVG)
# 2- DEG between clusters (DEG)
# 3- Others...

# This will be added to slot gene.list of scTypeEval object

get.HGV <- function(norm.mat,
                    ngenes = 500,
                    margin = 1L){
   # Here we will compute the variance across genes (rows) and select the HVGs
   gene_variance <- apply(norm.mat, margin, var)  # Variance per gene
   gene_mean <- apply(norm.mat, margin, mean)  # Mean expression per gene
   
   # Compute the coefficient of variation (CV) for each gene
   cv <- gene_variance / gene_mean
   
   # Rank genes by CV and select top n_hvg most variable genes
   hvg_genes <- order(cv, decreasing = TRUE)[1:ngenes]
   
   return(hvg_genes)
}

# get markers using scran findMarkers
get.DEG <- function(mat,
                    ident,
                    block = NULL,
                    ngenes = 500,
                    test.type = "t",
                    ncores = 1,
                    bparam = NULL,
                    progressbar = TRUE,
                    ...){
   
   param <- set_parallel_params(ncores = ncores,
                                bparam = bparam,
                                progressbar = progressbar)
   
   de <- scran::findMarkers(x = norm.mat,
                            group = ident,
                            pval.type = "all",
                            BPPARAM = param,
                            direction = "up",
                            ...)
   
   markers <- de |>
      dplyr::arrange(log2FC) |>
      head(ngenes)
      
      
   
   return(markers)
   
   
}


# gpFISH
get.gpsFISHMarkers <- function(sc_count,
                               ident,
                               rm_unannot = T,
                               total.ngenes = 500,
                               perCT.ngenes = NULL,
                               seed = 22,
                               ncores = 1){
   #https://htmlpreview.github.io/?https://github.com/kharchenkolab/gpsFISH/blob/main/doc/gene_panel_selection.html
   simulation_params <- gpsFISH::simulation_params
   
   # relative expression
   unique_cluster_label <- as.character(unique(ident))#unique cell type labels
   
   #cluster-wise relative expression
   relative_prop <-  list()
   relative_prop[["cluster.average"]] <- sapply(unique_cluster_label,
                                                gene_list=rownames(sc_count),
                                                gpsFISH::relative_freq,
                                                count_matrix=sc_count,             
                                                cell_cluster_conversion=ident)
   #individual cell level relative expression
   relative_prop[["cell.level"]] = t(t(sc_count)/colSums(sc_count))
   
   # Filter genes
   ave_count_ik_all_gene = sapply(unique_cluster_label,
                                  gene_list = rownames(sc_count),
                                  gpsFISH::average_count,
                                  count_matrix = sc_count,               
                                  cell_cluster_conversion = sc_cluster)
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
   
   gene2include.id <- which(rownames(sc_count) %in% gene2include.symbol)   #in the gene panel selection, genes will be encoded by their position in the gene list. Therefore, we need to get the position of curated marker genes here.
   gene2include.symbol <- rownames(sc_count)[gene2include.id]
   
   # Construct weighted penalty matrix based on cell type hierarchy
   cluster_distance <- gpsFISH::cluster_dis(count_table = sc_count,                       
                                            cell_cluster_conversion = sc_cluster,
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
   pop_size = 100 #population size: number of gene panels in a population
   panel_size = ngenes        #panel size: number of genes in a gene panel
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
