
dm_estimateMeanExpression <- function(counts, verbose = FALSE){
  
  # Calculate mean expression of genes 
  time_start <- Sys.time()
  if(verbose) message("* Calculating mean gene expression.. \n")
  
  inds <- 1:length(counts)
  
  mean_expression <- unlist(lapply(inds, function(g){ 
    mean(colSums(counts[[g]]), na.rm = TRUE)
  }))
  
  names(mean_expression) <- names(counts)
  
  time_end <- Sys.time()
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  return(mean_expression)
  
}


