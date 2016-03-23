#######################################################
#  group testing
#######################################################

#' @importFrom stats pchisq p.adjust

dmDS_test <- function(stats_full, stats_null, verbose = FALSE){
  
  ## calculate lr
  if(verbose) message("* Calculating likelihood ratio statistics.. \n")
  
  time_start <- Sys.time()
  
  lr <- 2*(rowSums(stats_full) - stats_null[, "lik"])
  
  nrgroups <- rowSums(!is.na(stats_full))
  
  df <- (nrgroups - 1) * stats_null[, "df"]
  
  df[nrgroups == 0] <- NA 
  lr[nrgroups == 0] <- NA 
  
  pvalue <- pchisq(lr, df = df , lower.tail = FALSE)
  
  adj_pvalue <- p.adjust(pvalue, method="BH")
  
  table <- data.frame(gene_id = rownames(stats_full), lr = lr, df = df, 
    pvalue = pvalue, adj_pvalue = adj_pvalue, stringsAsFactors = FALSE)
  
  # o <- order(table[, "pvalue"])
  # table <- table[o,]
  
  rownames(table) <- NULL
  
  time_end <- Sys.time()
  
  if(verbose) message("Took ", as.numeric(time_end - time_start), " seconds.\n")
  
  return(table)
  
}

