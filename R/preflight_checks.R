# function to purge sample and annotation names
purge_label <- function(label){
   label <- gsub(" |_|[+]|-", ".", label)
   label <- gsub(",", "", label)
   return(label)
}


check_ident <- function(scTypeEval = NULL,
                         ident,
                         verbose = TRUE){
   if(!is.null(scTypeEval)){
      if(is.null(ident)){
         ident <- scTypeEval@active_ident
         if(verbose){message("\nUsing default ident: ", ident, "\n")}
      }
      
      if(!ident %in% names(scTypeEval@metadata)){
         stop("Please provide a ident, i.e. a cell type or annotation to group cells included in metadata")
      }
      
      # retrieve ident and convert to factor
      ident_name <- ident
      ident <- scTypeEval@metadata[[ident]]
   }
   ident <- purge_label(ident)
   ident <- factor(ident)
   
   if(length(unique(ident)) < 2){
      stop("Less than 2 cell type detected in ident, Consistency metrics not possible")
   }
   
   return(ident)
}


check_sample <- function(scTypeEval = NULL,
                          sample,
                          verbose = TRUE){
   if(!is.null(scTypeEval)){
      if(!sample %in% names(scTypeEval@metadata)){
         stop("`sample` parameter not found in metadata colnames.\n")
      }
      # retrieve sample and convert to factor
      sample <- scTypeEval@metadata[[sample]]
   }
   sample <- purge_label(sample)
   sample <- factor(sample)
   
   if(length(levels(sample))<2){
      stop("For inter-sample comparison at least 2 samples is required, but more is recommended.\n")
   }
   
   if(length(levels(sample))<5){
      warning("Only ", length(levels(sample)), " samples detected.\n",
              "For inter-sample comparison 5 or more samples is recommended.\n")
   }
   
   return(sample)
}

check_blacklist <- function(scTypeEval, black_list = NULL, verbose = TRUE){
   
   if(is.null(black_list)){
      black_list <- scTypeEval@black_list
   }
   
   if(is.null(black_list)){
      if(verbose){message("Not using black gene list\n")}
   }
   return(black_list)
}

check_genelist <- function(scTypeEval, gene_list = NULL, verbose = TRUE){
   # set gene lists
   if(is.null(gene_list)){
      if(length(scTypeEval@gene_lists) == 0){
         stop("No gene list found.\n
              Add custom gene list or compute highly variable genes with add.HVG()\n")
      }
      
      gene_list <- names(scTypeEval@gene_lists)[1]
   } 
   
   if(!gene_list %in% names(scTypeEval@gene_lists)){
      stop("gene_list not included in scTypeEval object")
   }
   
   if(verbose){message("\nUsing ", gene_list, " gene list.\n")}
   gl <- scTypeEval@gene_lists[[gene_list]]
   
   if(length(gl)>5000){
      warning("Large number of genes selected for downstream analyses, this may increase computing time.\n")
   }
   
   return(gl)
}

check_dissimilarity_assays <- function(scTypeEval,
                                       slot = "all"){
   diss.assays <- names(scTypeEval@dissimilarity)
   if(length(diss.assays)<1){
      stop("No Dissimilarity slots found. Please run before `run_dissimilarity()`.\n")
   }
   
   if(length(slot)>1 || slot != "all"){
      diss.assays <- diss.assays[diss.assays %in% slot]
      if(length(diss.assays)<1){
         stop("No Dissimilarity slots found for ", slot, ". Please run before `run_dissimilarity()` or specific proper slot.\n")
      }
   }
   return(diss.assays)
}

check_dim_red_assays <- function(scTypeEval,
                                slot = "all"){
   red.assays <- names(scTypeEval@reductions)
   if(length(red.assays)<1){
      stop("No reduction slots found. Please run before `run_pca()` or `add_dim_reduction()` to add dimensional reductions slots.\n")
   }
   
   if(length(slot)>1 || slot != "all"){
      red.assays <- red.assays[red.assays %in% slot]
      if(length(red.assays)<1){
         stop("No reduction slots found for ", slot, ". Please run before `run_pca()` or `add_dim_reduction()` to add dimensional reductions slots.\n")
      }
   }
   return(red.assays)
}
