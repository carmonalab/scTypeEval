# function to purge sample and annotation names
purge_label <- function(label){
   label <- gsub(" |_|[+]|-", ".", label)
   label <- gsub(",", "", label)
   return(label)
}


.check_ident <- function(scTypeEval = NULL,
                         ident,
                         verbose = T){
   if(!is.null(scTypeEval)){
      if(is.null(ident)){
         ident <- scTypeEval@active.ident
         if(verbose){message("Using default ident: ", ident, "\n")}
      }
      
      if(!ident %in% names(scTypeEval@metadata)){
         stop("Please provide a ident, i.e. a cell type or annotation to group cells included in metadata")
      }
      
      # retrieve ident and convert to factor
      ident.name <- ident
      ident <- scTypeEval@metadata[[ident]]
   }
   ident <- purge_label(ident)
   ident <- factor(ident)
   
   if(length(unique(ident)) < 2){
      stop("Less than 2 cell type detected in ident, Consistency metrics not possible")
   }
   
   return(ident)
}


.check_sample <- function(scTypeEval = NULL,
                          sample,
                          verbose = T){
   if(!is.null(scTypeEval)){
      if(!sample %in% names(scTypeEval@metadata)){
         stop("`sample` parameter not found in metadata colnames.\n")
      }
      # retrieve sample and convert to factor
      sample <- scTypeEval@metadata[[sample]]
   }
   sample <- purge_label(sample)
   sample <- factor(sample)
   
   if(length(sample)<2){
      stop("For inter-sample comparison at least 2 samples is required, but more is recommended.\n")
   }
   
   if(length(sample)<5){
      warning("Only ", length(sample), " samples detected.\n",
              "For inter-sample comparison 5 or more samples is recommended.\n")
   }
   
   return(sample)
}

.check_blacklist <- function(scTypeEval, black.list, verbose = T){
   
   if(is.null(black.list)){
      black.list <- scTypeEval@black.list
   }
   
   if(is.null(black.list)){
      if(verbose){message("Not using black gene list\n")}
   }
   return(black.list)
}

.check_genelist <- function(scTypeEval, gene.list, verbose = TRUE){
   # set gene lists
   if(is.null(gene.list)){
      if(length(scTypeEval@gene.lists) == 0){
         stop("No gene list found to run consistency metrics.\n
              Add custom gene list or compute highly variable genes with add.HVG()\n")
      }
      
      gene.list <- names(scTypeEval@gene.lists)[1]
   } 
   
   if(gl %in% names(scTypeEval@gene.lists)){
      stop("Some gene list names not included in scTypeEval object")
   }
   
   if(verbose){message("Using ", gene.list, "gene list.\n")}
   gl <- scTypeEval@gene.lists[gene.list]
   
   if(length(gl)>5000){
      warning("Large number of genes selected for downstream analyses, this may increase computing time.\n")
   }
   
   return(gl)
}
