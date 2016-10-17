##############################################################################
# adjustements to profile likelihood
##############################################################################


dmSQTL_adjustmentOneGeneManyGroups <- function(g, counts, genotypes, pi, gamma0){  
  # g = 1; y = data@counts[[g]]; snps = data@genotypes[[g]]
  
  y <- counts[[g]]
  snps <- genotypes[[g]]
  adj <- rep(NA, nrow(snps))
  pig <- pi[[g]]
  
  for(i in 1:nrow(snps)){
    # i = 1
    
    if(all(is.na(pig[[i]][1, ]))) 
      next
    
    NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
    yg <- y[, !NAs]             
    group <- snps[i, !NAs]
    group <- factor(group)
    ngroups <- nlevels(group)
    lgroups <- levels(group)
    nlibs <- length(group)
    
    igroups <- lapply(lgroups, function(gr){which(group == gr)})
    names(igroups) <- lgroups
    
    adj[i] <- dm_adjustmentOneGeneManyGroups(y = yg, ngroups = ngroups, 
      lgroups = lgroups, igroups = igroups, gamma0 = gamma0, pi = pig[[i]]) 
    
  }
  return(adj)
  
}



dmSQTL_adjustmentCommon <- function(gamma0, counts, genotypes, pi, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  inds <- 1:length(counts)
  
  adj <- BiocParallel::bplapply(inds, dmSQTL_adjustmentOneGeneManyGroups, 
    counts = counts, genotypes = genotypes, pi = pi, gamma0 = gamma0, 
    BPPARAM = BPPARAM)
  
  
  adj <- unlist(adj)
  adj <- sum(adj, na.rm = TRUE)
  
  return(adj)
  
}














