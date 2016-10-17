##############################################################################
# Multiple group fitting 
# It returns a list of MatrixLists; unlistData = pi, metadata = stats
##############################################################################


dmSQTL_fitFull_gene <- function(g, counts, genotypes, 
  lgroups_g, ngroups_g, gamma0, prop_mode, prop_tol, verbose){
  # g = 9
  
  y <- counts[[g]]
  n_y <- nrow(y)
  snps <- genotypes[[g]]
  n_snps <- nrow(snps)
  names_snps <- rownames(snps)
  
  pi <- matrix(NA, nrow = n_y * n_snps, ncol = ngroups_g, 
    dimnames = list(rep(rownames(y), n_snps), lgroups_g))
  stats <- matrix(NA, n_snps, ngroups_g, 
    dimnames = list(names_snps, lgroups_g))
  
  for(i in 1:n_snps){          
    # i = 1
    
    NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
    yg <- y[, !NAs, drop = FALSE]             
    group <- snps[i, !NAs, drop = FALSE]
    group <- factor(group)
    ngroups <- nlevels(group)
    lgroups <- levels(group)
    nlibs <- length(group)
    
    igroups <- lapply(lgroups, function(gr){which(group == gr)})
    names(igroups) <- lgroups
    
    f <- dm_fitOneGeneManyGroups(y = yg, ngroups = ngroups, 
      lgroups = lgroups, igroups = igroups, 
      gamma0 = gamma0[[g]][i], prop_mode = prop_mode, 
      prop_tol = prop_tol, verbose = verbose)  
    
    ipi <- (i-1)*n_y + 1
    
    pi[ipi:(ipi+n_y-1), lgroups] <- f$pi
    stats[i, lgroups] <- f$stats
    
  }
  
  partitioning <- split(1:nrow(pi), factor(rep(1:n_snps, each = n_y)))
  names(partitioning) <- names_snps
  
  ff <- new("MatrixList", unlistData = pi, partitioning = partitioning, 
    metadata = stats)
  
  return(ff)
  
}

dmSQTL_fitNull_gene <- function(g, counts, genotypes, 
  gamma0, prop_mode, prop_tol, verbose){
  # g = 1; y = counts[[g]]; snps = genotypes[[g]]
  
  y = counts[[g]]
  n_y <- nrow(y)
  snps = genotypes[[g]]
  n_snps <- nrow(snps)
  names_snps <- rownames(snps)
  
  pi <- matrix(NA, nrow = n_y * n_snps, ncol = 1, 
    dimnames = list(rep(rownames(y), n_snps), "null"))
  stats <- matrix(NA, n_snps, 2, 
    dimnames = list(names_snps, c("lik", "df")))
  
  for(i in 1:n_snps){          
    # i = 1
    # print(i)
    
    NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
    yg <- y[, !NAs, drop = FALSE]
    
    f <- dm_fitOneGeneOneGroup(y = yg, gamma0 = gamma0[[g]][i], 
      prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
    
    ipi <- (i-1)*n_y + 1
    
    pi[ipi:(ipi+n_y-1), ] <- f$pi
    stats[i, ] <- f$stats
    
  }
  
  partitioning <- split(1:nrow(pi), factor(rep(1:n_snps, each = n_y)))
  names(partitioning) <- names_snps
  
  ff <- new("MatrixList", unlistData = pi, partitioning = partitioning,
    metadata = stats)
  
  return(ff)
  
}



dmSQTL_fitOneModel <- function(counts, genotypes, dispersion, model = "full", 
  prop_mode = "constrOptimG", prop_tol = 1e-12, verbose=FALSE, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  inds <- 1:length(counts)
  
  if(class(dispersion) == "numeric"){ 
    gamma0 <- relist(rep(dispersion, nrow(genotypes)), genotypes@partitioning)
  } else {
    gamma0 <- dispersion
  }
  
  lgroups_g <- c("0", "1", "2")
  ngroups_g <- 3
  
  switch(model, 
    
    full={
      
      if(verbose) message("* Fitting full model.. \n")
      
      time <- system.time(fff <- BiocParallel::bplapply(inds, dmSQTL_fitFull_gene, 
        counts = counts, genotypes = genotypes, lgroups_g = lgroups_g, 
        ngroups_g = ngroups_g, gamma0 = gamma0, prop_mode = prop_mode, 
        prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM))
      
      if(verbose) message("Took ", round(time["elapsed"], 2), " seconds.\n")
      names(fff) <- names(counts)  
      
      return(fff)
      
    }, 
    
    null={
      
      if(verbose) message("* Fitting null model.. \n")
      
      time <- system.time(fff <- BiocParallel::bplapply(inds, dmSQTL_fitNull_gene, 
        counts = counts, genotypes = genotypes, gamma0 = gamma0, 
        prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, 
        BPPARAM = BPPARAM))
      
      if(verbose) message("Took ", round(time["elapsed"], 2), " seconds.\n")
      names(fff) <- names(counts)
      
      return(fff)
      
    })
  
}

