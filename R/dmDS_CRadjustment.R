

dmDS_CRadjustmentManyGroups_gene <- function(g, counts, 
  ngroups, lgroups, igroups, 
  disp, prop, verbose){  
  # g = 1
  
  if(verbose >= 2)
    message(" Gene:", g)
  
  a <- dm_CRadjustmentManyGroups(y = counts[[g]], 
    ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
    disp = disp[g], prop = prop[[g]])
  
  return(a)
  
}



dmDS_CRadjustment <- function(counts, fit, design, dispersion,
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  if(verbose) message("* Calculating Cox-Reid adjustment.. \n")
  
  inds <-  1:length(counts)
  
  # Prepare dispersion
  if(length(dispersion) == 1){
    disp <- rep(dispersion, length(inds))
  } else {
    disp <- dispersion
  }
  
  # If the design is equivalent to a oneway layout, use a shortcut algorithm
  groups <- edgeR::designAsFactor(design)
  
  if(nlevels(groups) == ncol(design)){
    
    groups <- factor(groups, labels = paste0("gr", levels(groups)))
    ngroups <- nlevels(groups)
    lgroups <- levels(groups)
    igroups <- lapply(lgroups, function(gr){which(groups == gr)})
    names(igroups) <- lgroups
    
    figroups <- unlist(lapply(igroups, function(x){x[1]}))
    
    prop <- fit[, figroups]
    
    time <- system.time(aa <- BiocParallel::bplapply(inds, 
      dmDS_CRadjustmentManyGroups_gene, counts = counts, 
      ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
      disp = disp, prop = prop, 
      verbose = verbose, BPPARAM = BPPARAM))
    
    names(aa) <- names(counts)
    
    adj <- unlist(aa) 
    
  }else{
    stop("Currently, regression framework is not implemented here!")
  }
  
  if(verbose >= 2) message("\n")
  if(verbose) message("Took ", round(time["elapsed"], 4), " seconds.\n")
  
  # adj is a vector
  return(adj)
  
  
}












