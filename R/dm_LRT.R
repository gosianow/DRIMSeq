#' @importFrom stats pchisq p.adjust


dm_LRT <- function(lik_full, lik_null, df, verbose = FALSE){
  
  if(verbose) 
    message("* Calculating likelihood ratio statistics.. \n")
  
  time_start <- Sys.time()
  
  lr <- 2*(lik_full - lik_null)
  
  pvalue <- pchisq(lr, df = df , lower.tail = FALSE)
  
  adj_pvalue <- p.adjust(pvalue, method="BH")
  
  table <- data.frame(lr = lr, df = df, 
    pvalue = pvalue, adj_pvalue = adj_pvalue, 
    stringsAsFactors = FALSE)

  rownames(table) <- names(lik_full)
  
  time_end <- Sys.time()
  
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  return(table)
  
}



