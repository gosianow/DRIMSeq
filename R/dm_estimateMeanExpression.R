##############################################################################
# calculate common dispersion 
##############################################################################

dm_estimateMeanExpression <- function(counts, verbose = FALSE){
  
  ### calculate mean expression of genes 
  if(verbose) message("* Calculating mean gene expression.. \n")
  
  inds <- 1:length(counts)
  
  time <- system.time(mean_expression <- unlist(lapply(inds, function(g){ 
    
    mean(colSums(counts[[g]]), na.rm = TRUE)
    
  })))
  
  names(mean_expression) <- names(counts)
  
  if(verbose) message("Took ", round(time["elapsed"], 2), " seconds.\n")
  
  return(mean_expression)
  
}


