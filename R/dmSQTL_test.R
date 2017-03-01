#######################################################
#  group testing
#######################################################

#' @importFrom stats pchisq p.adjust



dmSQTL_test_per_gene <- function(g, fit_full, fit_null, gene_list){
  # g = 662
  
  lr <- 2*(rowSums(fit_full[[g]]@metadata, na.rm = TRUE) - 
      fit_null[[g]]@metadata[, "lik"])
  
  nrgroups <- rowSums(!is.na(fit_full[[g]]@metadata))
  
  ### negative when NAs in all groups in lik
  df <- (nrgroups - 1)*fit_null[[g]]@metadata[, "df"] 
  
  df[nrgroups == 0] <- NA 
  lr[nrgroups == 0] <- NA 
  
  pvalue <- pchisq(lr, df = df , lower.tail = FALSE)
  
  tt <- data.frame(gene_id = gene_list[g], 
    snp_id = rownames(fit_full[[g]]@metadata), lr = lr, df = df, 
    pvalue = pvalue, stringsAsFactors = FALSE)
  
}




dmSQTL_test <- function(fit_full, fit_null, return_list = FALSE, 
  verbose = FALSE, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  ## calculate lr
  if(verbose) message("* Calculating likelihood ratio statistics.. \n")
  time_start <- Sys.time()
  
  inds <- 1:length(fit_full)
  gene_list <- names(fit_full)
  
  table_list <- BiocParallel::bplapply(inds, dmSQTL_test_per_gene, 
    fit_full = fit_full, fit_null = fit_null, gene_list = gene_list, 
    BPPARAM = BPPARAM)
  
  if(return_list)
    return(table_list)
  
  table <- do.call(rbind, table_list)
  
  adj_pvalue <- p.adjust(table[, "pvalue"], method="BH")
  
  table$adj_pvalue <- adj_pvalue
  
  # o <- order(table[, "pvalue"])  
  # table <- table[o,]
  
  rownames(table) <- NULL
  
  time_end <- Sys.time()
  if(verbose) message("Took ", as.numeric(time_end - time_start), " seconds.\n")
  
  return(table)
  
}


