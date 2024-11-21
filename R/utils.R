# S4 class object
methods::setClass("scTypeEval",
                  slots = c(
                     counts = "dgCMatrix", # raw counts
                     metadata = "data.frame", # data.frame with metadata
                     norm.param = "list", # normalization paramters to normalize a subset of genes on the fly
                     consistency = "list", # actual consistency results assays
                     gene.lists = "list", # list of either HGV, markers
                     black.list = "character", # list of genes in black list
                     active.ident = "factor", # default grouping variable
                     reductions = "list", # list of dim reductions
                     misc = "list", # miscellaneous
                     commands = "list", # commands performed
                     version = "character" # package version
                  )
)

# Define the initialize method to set default for reductions and consistency
methods::setMethod("initialize", "scTypeEval", function(.Object, ...) {
   args <- list(...)
   
   # Set default for consistency if not provided
   if (is.null(args$consistency)) {
      args$consistency <- list(intraSample = NULL, interSample = NULL)
   }
   
   # Set default for reductions if not provided
   if (is.null(args$reductions)) {
      args$reductions <- list(intraSample = NULL, interSample = NULL)
   }
   
   # Set default for metadata if not provided
   if (is.null(args$metadata)) {
      args$metadata <- S4Vectors::DataFrame()
   }
   
   # Pass the updated arguments to the default initialize method
   .Object <- callNextMethod(.Object, ..., 
                             consistency = args$consistency, 
                             reductions = args$reductions, 
                             metadata = args$metadata)
   
   validObject(.Object)  # Validate the object
   .Object
})

# define consistency assay object
methods::setClass("ConsistencyAssay",
                  slots = c(
                     metrics = "data.frame",
                     method = "character",
                     gene.list = "character",
                     black.list = "character",
                     ident = "factor"
                  )
)

methods::setClass("DimReduc",
                  slots = c(
                     cell.embeddings = 'matrix',
                     feature.loadings = 'matrix',
                     feature.loadings.projected = 'matrix',
                     gene.list = "character",
                     black.list = "character",
                     key = "character" # type of reduction, PCA, UMAP...
                  )
)





set_parallel_params <- function(ncores = NULL,
                                bparam = NULL,
                                progressbar = TRUE)
{
   if (is.null(ncores)) {
      ncores <- 1
   }
   
   if (ncores > parallel::detectCores()) {
      ncores <- parallel::detectCores()
      message("Using more cores available in this computer, reducing number of cores to ",
              ncores)
   }
   
   # set parallelization parameters
   if (is.null(bparam)) {
      if (ncores > 1) {
         if(.Platform$OS.type == "windows"){
            param <- BiocParallel::SnowParam(workers=ncores,
                                             progressbar = progressbar)
         } else {
            param <- BiocParallel::MulticoreParam(workers =  ncores,
                                                  progressbar = progressbar)
         }
      } else {
         param <- BiocParallel::SerialParam(progressbar = progressbar)
      }
   } else {
      param <- bparam
   }
   return(param)
}

require(ggplot2)
single.plot.consistency <- function(tbl,
                                    metric = "modularity",
                                    only.mean = F){
   
   tit <- unique(tbl$pseudomodality)
   
   pl <- tbl %>%
      ggplot() +
      geom_line(aes(true_ncells, mean,
                    color = Shuffling_rate,
                    linetype = Shuffling_method)) +
      scale_linetype_manual(values = c(1,6,3)) +
      labs(y = paste(metric, "score") ,
           x = "granularity (n cells)",
           title = tit,
           color = "Agreement after\nlabel shuffling") +
      ggprism::theme_prism() +
      theme(legend.title = element_text(size = 12,
                                        face = "bold"))
   
   if(!only.mean){
      pl <- pl +
         geom_jitter(aes(true_ncells, .data[[metric]],
                         color = Shuffling_rate),
                     alpha = 0.2,
                     width = 0.75) +
         facet_wrap(~ Shuffling_rate + Shuffling_method,
                    ncol = 2)
      
   }
   
   
   return(pl)
   
}


single.plot.distrCon <- function(tbl,
                                 metric = "modularity",
                                 only.mean = F){
   
   
   
   if(only.mean){
      tit <- unique(tbl$pseudomodality)
      # do average
      tbl <- tbl %>%
         group_by(instance) %>%
         summarize(mean = mean(.data[[metric]][celltype != "unknown"]),
                   sd = sd(.data[[metric]][celltype != "unknown"])
         ) %>%
         ungroup()
      
      tbl[["val"]] <- tbl[["mean"]]
      tbl <- list(tbl)
      
   } else {
      tbl[["val"]] <- tbl[[metric]]
      tbl <- split(tbl, tbl$celltype)
      tit <- NULL
   }
   
   pls <- lapply(tbl, function(t){
      if(is.null(tit)){
         tit <- unique(t$celltype)
      }
      
      d1 <- t %>% filter(instance != "original")
      d2 <- t %>% filter(instance == "original")
      
      p_zval <- p_val_zscore(d2[["val"]], d1[["val"]])
      
      tit2 <- paste(tit, "\nZ score:",
                    round(p_zval$z_score, 1),
                    "/ p-val:",
                    format(p_zval$p_val, scientifc = T, digits = 1)
      )
      
      p <- ggplot() +
         geom_density(data = d1,
                      aes(x = val),
                      color = "grey15",
                      fill = "grey70",
                      alpha = 0.6) +
         geom_vline(xintercept = d2[["val"]],
                    linetype = 2,
                    color = "red") +
         xlim(c(-0.35, 1)) +
         labs(title = tit2,
              x = paste0(metric, " score")
         ) +
         ggprism::theme_prism()
      return(p)
   })
   
   
   if(length(pls) == 1){
      pls <- pls[[1]]
   } else {
      pls <- ggpubr::ggarrange(plotlist = pls)
   }
   
   return(pls)
}


fun.sel <- function(x, version){
   x %>%
      group_by(pseudomodality) %>%
      summarize(silhouette = mean(silhouette),
                n_cells = first(n_cells)) %>%
      ungroup() %>%
      mutate(version = version)
}


p_val_zscore <- function(obs,
                         random_values) {
   
   mean <- mean(random_values)
   sd <- sd(random_values)
   z_score <- (obs - mean) / sd
   p_val <- 1 - stats::pnorm(z_score)
   
   ret <- list(z_score = z_score,
               p_val = p_val)
   return(ret)
}


get_intra_SvB <- function(object,
                          colLabel,
                          celltype,
                          min.cells = 10,
                          seed = 22,
                          assay = "RNA",
                          ncores = 1) {
   
   # set paralelization
   param <- set_parallel_params(ncores = ncores, progressbar = F)
   
   ct.mat <- object[, object@meta.data[[colLabel]] == celltype][[assay]]$counts
   colnames(ct.mat) <- gsub("_", "-", colnames(ct.mat)) %>%
      paste(., celltype, sep = "_")
   
   ncells <- ncol(ct.mat)
   
   # adjust if their proportions is very discrepant
   seu.psblk <- object[, object@meta.data[[colLabel]] != celltype]
   
   # find variable featurs
   seu.psblk <- FindVariableFeatures(seu.psblk,
                                     selection.method = "vst",
                                     nfeatures = 1000,
                                     verbose = F)
   
   # produce sketching and merge recursively
   ps.mat <- bplapply(1:ncells,
                      BPPARAM = param,
                      function(n){
                         # sketching
                         
                         myseed <- seed + n
                         suppressMessages(
                            {
                               ds <- Seurat::SketchData(seu.psblk,
                                                        assay = assay,
                                                        ncells = min.cells,
                                                        seed = myseed,
                                                        verbose = F)
                               ds.mat <- ds[["sketch"]]$counts
                               ds.mat <- Matrix::rowSums(ds.mat)
                            })
                         
                         return(ds.mat)
                      })
   
   # join Sketch pseudobulk
   ps.mat <- do.call(cbind, ps.mat)
   
   colnames(ps.mat) <- paste0(paste0("r", 1:ncells), "_psblk")
   
   # join with actual celltypes
   fi <- cbind(ct.mat, ps.mat)
   
   return(fi)
   
   
}


clr_function <- function(x){
   #https://github.com/satijalab/seurat/blob/1549dcb3075eaeac01c925c4b4bb73c73450fc50/R/preprocessing.R#L4418
   return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
}
