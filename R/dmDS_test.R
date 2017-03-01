#######################################################
#  group testing
#######################################################

#' @importFrom stats pchisq p.adjust

dmDS_test <- function(lik_full, lik_null, df, verbose = FALSE){
  
  ## calculate lr
  if(verbose) message("* Calculating likelihood ratio statistics.. \n")
  
  time_start <- Sys.time()
  
  lr <- as.numeric(2*(rowSums(lik_full) - lik_null))
  
  nrgroups <- rowSums(!is.na(lik_full))
  
  df <- (nrgroups - 1) * df
  
  df[nrgroups == 0] <- NA 
  lr[nrgroups == 0] <- NA 
  
  pvalue <- pchisq(lr, df = df , lower.tail = FALSE)
  
  adj_pvalue <- p.adjust(pvalue, method="BH")
  
  table <- data.frame(lr = lr, df = df, 
    pvalue = pvalue, adj_pvalue = adj_pvalue, 
    stringsAsFactors = FALSE)

  rownames(table) <- rownames(lik_full)
  
  time_end <- Sys.time()
  
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  return(table)
  
}

