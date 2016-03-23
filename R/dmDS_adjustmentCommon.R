##############################################################################
# adjustements to profile likelihood for common dispersion -> sum
##############################################################################


dmDS_adjustmentOneGeneManyGroups <- function(g, counts, 
  ngroups, lgroups, igroups, gamma0, pi){  
  # g = 1
  
  a <- dm_adjustmentOneGeneManyGroups(y = counts[[g]], ngroups = ngroups, 
    lgroups = lgroups, igroups = igroups, gamma0 = gamma0, pi = pi[[g]]) 
  
  return(a)
  
}



dmDS_adjustmentCommon <- function(gamma0, counts, samples, pi, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  inds <- 1:length(counts)
  
  group <- samples$group
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  
  igroups <- lapply(lgroups, function(gr){which(group == gr)})
  names(igroups) <- lgroups
  
  adj <- BiocParallel::bplapply(inds, dmDS_adjustmentOneGeneManyGroups, 
    counts = counts, ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
    gamma0 = gamma0, pi = pi, BPPARAM = BPPARAM)
  
  adj <- unlist(adj)
  adj <- sum(adj, na.rm = TRUE) 
  
  return(adj)
  
}

