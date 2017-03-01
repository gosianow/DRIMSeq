#######################################################
#  group testing
#######################################################


#' Perform the two stage test
#' 
#' @param pvalue_gene data frame with pvalue and gene_id
#' @param pvalue_feature data frame with pvalue, gene_id and feature_id
#' @param FDR numeric cutoff for the FDR
#' @param verbose logical 
#' @return Returns a data frame with adjusted feature level p-values
#' @importFrom stats pchisq p.adjust
#' 
dm_twoStageTest <- function(pvalue_gene, pvalue_feature, FDR = 0.05, 
  verbose = FALSE){
  
  if(verbose) message("* Perform the two stage test.. \n")
  
  time_start <- Sys.time()
  
  pvalue_gene[, "adj_pvalue"] <- p.adjust(pvalue_gene[, "pvalue"], method="BH")
  
  genes2keep <- pvalue_gene[
    pvalue_gene[, "adj_pvalue"] < FDR & !is.na(pvalue_gene[, "adj_pvalue"]), 
    "gene_id", drop = FALSE]
  
  pvalue_feature_split <- split(pvalue_feature[, "pvalue"], 
    factor(pvalue_feature[, "gene_id"], 
    levels = pvalue_gene[, "gene_id"]))
  
  pvalue_two_stages <- lapply(names(pvalue_feature_split), function(i){
    
    x <- pvalue_feature_split[[i]]
    
    if(i %in% genes2keep){
      out <- p.adjust(x, method="bonferroni")
    }else{
      out <- rep(1, length(x))
    }
    
    return(out)
  })
  
  pvalue_feature[, "adj_pvalue"] <- unlist(pvalue_two_stages)
  
  time_end <- Sys.time()
  
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  return(pvalue_feature[, c("gene_id", "feature_id", "adj_pvalue")])
  
}





















