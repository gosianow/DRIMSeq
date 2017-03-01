##############################################################################
# Multiple group fitting 
# It returns a list of MatrixLists; unlistData = prop, metadata = stats
##############################################################################


dmSQTL_fitFull_gene <- function(g, counts, genotypes, 
  lgroups_g, ngroups_g, disp, prop_mode, prop_tol, verbose){
  # g = 9
  
  y <- counts[[g]]
  n_y <- nrow(y)
  snps <- genotypes[[g]]
  n_snps <- nrow(snps)
  names_snps <- rownames(snps)
  
  prop <- matrix(NA, nrow = n_y * n_snps, ncol = ngroups_g, 
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
      disp = disp[[g]][i], prop_mode = prop_mode, 
      prop_tol = prop_tol, verbose = verbose)  
    
    iprop <- (i-1)*n_y + 1
    
    prop[iprop:(iprop+n_y-1), lgroups] <- f$prop
    stats[i, lgroups] <- f$stats
    
  }
  
  partitioning <- split(1:nrow(prop), factor(rep(1:n_snps, each = n_y)))
  names(partitioning) <- names_snps
  
  ff <- new("MatrixList", unlistData = prop, partitioning = partitioning, 
    metadata = stats)
  
  return(ff)
  
}

dmSQTL_fitNull_gene <- function(g, counts, genotypes, 
  disp, prop_mode, prop_tol, verbose){
  # g = 1; y = counts[[g]]; snps = genotypes[[g]]
  
  y = counts[[g]]
  n_y <- nrow(y)
  snps = genotypes[[g]]
  n_snps <- nrow(snps)
  names_snps <- rownames(snps)
  
  prop <- matrix(NA, nrow = n_y * n_snps, ncol = 1, 
    dimnames = list(rep(rownames(y), n_snps), "null"))
  stats <- matrix(NA, n_snps, 2, 
    dimnames = list(names_snps, c("lik", "df")))
  
  for(i in 1:n_snps){          
    # i = 1
    # print(i)
    
    NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
    yg <- y[, !NAs, drop = FALSE]
    
    f <- dm_fitOneGeneOneGroup(y = yg, disp = disp[[g]][i], 
      prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
    
    iprop <- (i-1)*n_y + 1
    
    prop[iprop:(iprop+n_y-1), ] <- f$prop
    stats[i, ] <- f$stats
    
  }
  
  partitioning <- split(1:nrow(prop), factor(rep(1:n_snps, each = n_y)))
  names(partitioning) <- names_snps
  
  ff <- new("MatrixList", unlistData = prop, partitioning = partitioning,
    metadata = stats)
  
  return(ff)
  
}



dmSQTL_fitOneModel <- function(counts, genotypes, dispersion, model = "full", 
  prop_mode = "constrOptim", prop_tol = 1e-12, verbose=FALSE, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  inds <- 1:length(counts)
  
  if(class(dispersion) == "numeric"){ 
    disp <- relist(rep(dispersion, nrow(genotypes)), genotypes@partitioning)
  } else {
    disp <- dispersion
  }
  
  lgroups_g <- c("0", "1", "2")
  ngroups_g <- 3
  
  switch(model, 
    
    full={
      
      if(verbose) message("* Fitting full model.. \n")
      
      time <- system.time(fff <- BiocParallel::bplapply(inds, dmSQTL_fitFull_gene, 
        counts = counts, genotypes = genotypes, lgroups_g = lgroups_g, 
        ngroups_g = ngroups_g, disp = disp, prop_mode = prop_mode, 
        prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM))
      
      if(verbose) message("Took ", round(time["elapsed"], 2), " seconds.\n")
      names(fff) <- names(counts)  
      
      return(fff)
      
    }, 
    
    null={
      
      if(verbose) message("* Fitting null model.. \n")
      
      time <- system.time(fff <- BiocParallel::bplapply(inds, dmSQTL_fitNull_gene, 
        counts = counts, genotypes = genotypes, disp = disp, 
        prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, 
        BPPARAM = BPPARAM))
      
      if(verbose) message("Took ", round(time["elapsed"], 2), " seconds.\n")
      names(fff) <- names(counts)
      
      return(fff)
      
    })
  
}

