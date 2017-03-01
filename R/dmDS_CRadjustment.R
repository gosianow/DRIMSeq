

dmDS_CRadjustmentManyGroups_gene <- function(g, counts, 
  ngroups, lgroups, igroups, prec, prop, verbose){  

  if(verbose >= 2) message(" Gene:", g)
  
  a <- dm_CRadjustmentManyGroups(y = counts[[g]], 
    ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
    prec = prec[g], prop = prop[[g]])
  
  return(a)
  
}


dmDS_CRadjustmentRegression_gene <- function(g, counts, 
  design, prec, fit, verbose){  
	# g = 3
  if(verbose >= 2) message(" Gene:", g)
  
  a <- dm_CRadjustmentRegression(y = counts[[g]], x = design, 
    prec = prec[g], prop = fit[[g]])
  
  return(a)
  
}



dmDS_CRadjustment <- function(counts, fit, design, precision,
  one_way = TRUE,
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  time_start <- Sys.time()
  if(verbose) message("* Calculating Cox-Reid adjustment.. \n")
  
  inds <-  1:length(counts)
  
  # Prepare precision
  if(length(precision) == 1){
    prec <- rep(precision, length(inds))
  } else {
    prec <- precision
  }
  
  # If the design is equivalent to a oneway layout, use a shortcut algorithm
  groups <- edgeR::designAsFactor(design)
  
  if(nlevels(groups) == ncol(design) && one_way){
    
    groups <- factor(groups, labels = paste0("gr", levels(groups)))
    ngroups <- nlevels(groups)
    lgroups <- levels(groups)
    igroups <- lapply(lgroups, function(gr){which(groups == gr)})
    names(igroups) <- lgroups
    # Get the column number of a first occurance of a group level
    figroups <- unlist(lapply(igroups, function(x){x[1]}))
    
    prop <- fit[, figroups]
    
    a <- BiocParallel::bplapply(inds, 
      dmDS_CRadjustmentManyGroups_gene, counts = counts, 
      ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
      prec = prec, prop = prop, 
      verbose = verbose, BPPARAM = BPPARAM)
    
    names(a) <- names(counts)
    adj <- unlist(a) 
    
  }else{
    
    a <- BiocParallel::bplapply(inds, 
      dmDS_CRadjustmentRegression_gene, counts = counts, 
      design = design, prec = prec, fit = fit, 
      verbose = verbose, BPPARAM = BPPARAM)
    
    names(a) <- names(counts)
    adj <- unlist(a) 

  }
  
  time_end <- Sys.time()
  if(verbose >= 2) message("\n")
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  # adj is a vector of length G
  return(adj)

}












