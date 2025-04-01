############################################################################################################
# scAnnot: an automated annotation package for single-cell analysis
# Version: V1.0.0
# Creator: Lu Pan
# Last Update Date: 2025-04-01
############################################################################################################
#' @import fgsea
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Automated annotation package for single-cell analysis
#'
#' This function will run automated annotation for single-cell data.
#'
#' @param project_name Project name. 'scAnnot_Project' by default.
#' @param input_dir Input directory. Default directory
#' is the current working directory.
#' @param output_dir Output directory. Default directory is the
#' current working directory.
#' @param ref_degs_path Path to reference celltype-specific
#' differentially expressed genes (DEGs). Input format should be CSV or txt.
#' @param ref_degs Data frame of the reference celltype-specific degs.
#' @param ref_degs_format Data format of the reference celltype-specific degs.
#' Default is 'general'. Takes input: 'general', 'seurat', 'scanpy'.
#' 'general' format should include data columns: 'group','names','logfoldchanges',
#' 'pvals_adj','pvals', for which 'group' refers to the column for celltype names,
#' and 'names' refers to the column for gene names.
#' @param ref_padj_cutoff P-adjusted threshold. Default is 0.05.
#' @param ref_logfc_cutoff Log(fold-change) threshold. Default is 0.25.
#' @param ref_topn Top n genes from each reference celltype DEGs to be used
#' for annotation. Default is 30.
#' @param data_degs_path Path to cluster-specific DEGs of the input data to be annotated.
#' Input format should be CSV or txt.
#' @param data_degs Data frame of the cluster-specific degs of the input data.
#' @param data_degs_format Data format of the cluster-specific degs of the input data.
#' Default is 'general'. Takes input: 'general', 'seurat', 'scanpy'.
#' 'general' format should include data columns: 'group','names','logfoldchanges',
#' 'pvals_adj','pvals', for which 'group' refers to the column for celltype names,
#' and 'names' refers to the column for gene names.
#' @param data_padj_cutoff P-adjusted threshold of the cluster-specific degs
#' of the input data. Default is 0.05.
#' @param result_padj_cutoff P-adjusted threshold of the final result.
#' Default is 0.1.
#' @param filter Whether to filter and keep only the top predicted celltype for each cluster.
#' Default is set to TRUE.
#' @export

scAnnot <- function(project_name = "scAnnot_Project",
                      input_dir = "./",
                      output_dir = "./",
                      ref_degs_path = NULL, # path to reference celltype-specific degs
                      ref_degs = NULL, # data frame of the reference celltype-specific degs
                      ref_degs_format = "general", # general, seurat, scanpy
                      ref_padj_cutoff = 0.05,
                      ref_logfc_cutoff = 0.25,
                      ref_topn = 30,
                      data_degs_path = NULL, # path to data celltype-specific degs
                      data_degs = NULL, # data frame of the data celltype-specific degs
                      data_degs_format = "general", # general, seurat, scanpy
                      data_padj_cutoff = 0.05,
                      result_padj_cutoff = 0.1,
                      filter = TRUE){
  print("Initialising pipeline environment..")
  if(is.null(ref_degs_path)){
    if(is.null(ref_degs)){
      stop("Reference DEGs are not provided!")
    }
  }else{
    if(!file.exists(ref_degs_path)){
      stop("Path to reference DEGs is incorrect! Check your file path.")
    }
  }
  
  if(file.exists(ref_degs_path)){
    ref_degs <- read.csv(ref_degs_path, header = T)
    if(ncol(ref_degs) < 2){
      ref_degs <- read.table(ref_degs_path, sep = "\t", header = T)
    }
  }
  
  ref_padj <- NULL
  ref_logfc <- NULL
  ref_group <- NULL
  ref_genes <- NULL
  
  if(toupper(ref_degs_format) %in% c("SCANPY","GENERAL")){
    ref_padj <- "pvals_adj"
    ref_logfc <- "logfoldchanges"
    ref_group <- "group"
    ref_genes <- "names"
  }else if(toupper(ref_degs_format) == "SEURAT"){
    ref_padj <- "p_val_adj"
    ref_logfc <- "avg_log2FC"
    ref_group <- "cluster"
    ref_genes <- "gene"
  }

  ref_degs <- ref_degs[which(ref_degs[,ref_padj] < ref_padj_cutoff),]
  ref_degs <- split(ref_degs, ref_degs[,ref_group])
  
  ref <- NULL
  for(i in 1:length(ref_degs)){
    print(i)
    cdegs <- NULL
    cdegs <- ref_degs[[i]]
    cname <- NULL
    cname <- names(ref_degs)[i]
    if(nrow(cdegs) == 0){
      stop(paste("Provided reference ",cname," has no DEGs based on the selected p-adjusted cut-off. Please re-enter a higher p-adjusted cut-off threshold or submit other reference DEGs.", sep = ""))
    }
    ref[[i]] <- cdegs[1:ifelse(nrow(cdegs) > ref_topn, ref_topn, nrow(cdegs)),ref_genes]
    names(ref)[i] <- cname
  }

  if(is.null(data_degs_path)){
    if(is.null(data_degs)){
      stop("Data DEGs are not provided!")
    }
  }else{
    if(!file.exists(data_degs_path)){
      stop("Path to data DEGs is incorrect! Check your file path.")
    }
  }
  
  if(file.exists(data_degs_path)){
    data_degs <- read.csv(data_degs_path, header = T)
    if(ncol(data_degs) < 2){
      data_degs <- read.table(data_degs_path, sep = "\t", header = T)
    }
  }
  
  data_padj <- NULL
  data_logfc <- NULL
  data_group <- NULL
  data_genes <- NULL
  
  if(toupper(data_degs_format) %in% c("SCANPY","GENERAL")){
    data_padj <- "pvals_adj"
    data_logfc <- "logfoldchanges"
    data_group <- "group"
    data_genes <- "names"
  }else if(toupper(data_degs_format) == "SEURAT"){
    data_padj <- "p_val_adj"
    data_logfc <- "avg_log2FC"
    data_group <- "cluster"
    data_genes <- "gene"
  }
  
  
  data_degs <- data_degs[which(data_degs[,data_padj] < data_padj_cutoff),]
  data_degs <- split(data_degs, data_degs[,data_group])

  print("Running annotation based on provided reference..")
  
  final <- NULL
  for(i in 1:length(data_degs)){
    print(i)
    subcin <- NULL
    subcin <- data_degs[[i]]
    if(nrow(subcin) == 0){
      stop(paste("Provided data ",names(data_degs)[i]," has no DEGs based on the selected p-adjusted cut-off. Please re-enter a higher p-adjusted cut-off threshold or submit other data DEGs.", sep = ""))
    }
    subcin <- subcin[order(subcin[,data_logfc], decreasing = T),]
    crank <- NULL
    crank <- subcin[,data_logfc]
    names(crank) <- subcin[,data_genes]
    cresult <- NULL
    cresult <- fgseaMultilevel(pathways = ref, stats = crank, nPermSimple = 100, scoreType = ifelse(length(which(subcin[,data_logfc] < 0)) > 0, "std", "pos"))
    final <- rbind(final, data.frame(cluster = names(data_degs)[i], data.frame(cresult)))
  }
  
  final <- final[which(final$NES > 0),]
  final <- split(final, final$cluster)
  final <- lapply(final, function(x){
    x <- x[order(x$NES, decreasing = T),]
    return(x)
  })
  final <- do.call(rbind.data.frame, final)
  final <- final[which(final$padj < result_padj_cutoff),]
  row.names(final) <- NULL
  
  if(filter == TRUE){
    final <- split(final, final$cluster)
    final <- lapply(final, function(x){
      x <- x[order(x$NES, decreasing = T),]
      x <- x[1:ifelse(nrow(x) > 1, 1, nrow(x)),]
      return(x)
    })
    final <- do.call(rbind.data.frame, final)
    row.names(final) <- NULL
  }
  
  if(length(unique(final$cluster)) < length(names(data_degs))){
    final <- rbind(final, data.frame(cluster = names(data_degs)[which(!names(data_degs) %in% unique(final$cluster))], pathway = "Unknown", pval = NA, padj = NA, log2err = NA, ES = NA, NES = NA, size = NA, leadingEdge = NA))
  }
  
  colnames(final)[which(colnames(final) == "pathway")] <- "Celltype"
  row.names(final) <- NULL
  saveRDS(final, paste(output_dir,"/",project_name,"_scAnnot_Annotation_Results_Padj",result_padj_cutoff,"_.txt",sep = ""))
  
  final <- data.frame(apply(final,2,as.character))
  write.table(final, paste(output_dir,"/",project_name,"_scAnnot_Annotation_Results_Padj",result_padj_cutoff,".txt",sep = ""), row.names = F, sep = "\t", quote = F)

  print("Done!")
  return(final)
}
