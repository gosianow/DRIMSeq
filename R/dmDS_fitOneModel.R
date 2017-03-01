##############################################################################
# multiple group ffting 
##############################################################################

dmDS_fitFull_gene <- function(g, counts, 
  ngroups, lgroups, igroups, gamma0, prop_mode, prop_tol, verbose){  
  # g = 1
  if(verbose >= 2)
    message(" Gene:", g)
  
  f <- dm_fitOneGeneManyGroups(y = counts[[g]], ngroups = ngroups, 
    lgroups = lgroups, igroups = igroups, gamma0 = gamma0[g], 
    prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
  
  return(f)
  
}


dmDS_fitNull_gene <- function(g, counts, 
  gamma0, prop_mode, prop_tol, verbose){  
  # g = 1
  # message("Gene:", g)
  
  f <- dm_fitOneGeneOneGroup(y = counts[[g]], gamma0 = gamma0[g], 
    prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
  
  return(f)
  
}



dmDS_fitOneModel <- function(counts, samples, dispersion, model = "full",
  prop_mode = "constrOptim2G", prop_tol = 1e-12, verbose = FALSE, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  inds <-  1:length(counts)
  
  if(length(dispersion) == 1){
    gamma0 <- rep(dispersion, length(inds))
  } else {
    gamma0 <- dispersion
  }
  
  switch(model, 
    
    full = {
      
      if(verbose) message("* Fitting full model.. \n")
      
      group <- samples$group
      ngroups <- nlevels(group)
      lgroups <- levels(group)
      
      igroups <- lapply(lgroups, function(gr){which(group == gr)})
      names(igroups) <- lgroups
      
      time <- system.time(ff <- BiocParallel::bplapply(inds, dmDS_fitFull_gene, 
        counts = counts, ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
        gamma0 = gamma0, prop_mode = prop_mode, prop_tol = prop_tol, 
        verbose = verbose, BPPARAM = BPPARAM))
      
      names(ff) <- names(counts)
      
      fit <- MatrixList(lapply(ff, function(f) f[["pi"]]))
      
      lik <- do.call(rbind, lapply(ff, function(f) f[["lik"]])) 
      rownames(lik) <- names(counts)
      
      df <- rep(NA, nrow(lik))
      names(df) <- names(counts)
      
      if(verbose >= 2) message("\n")
      if(verbose) message("Took ", round(time["elapsed"], 2), " seconds.\n")
      
      return(list(fit = fit, lik = lik, df = df))
      
    },
    
    null = {
      
      if(verbose) message("* Fitting null model.. \n")
      
      time <- system.time(ff <- BiocParallel::bplapply(inds, dmDS_fitNull_gene, 
        counts = counts, gamma0 = gamma0, prop_mode = prop_mode, 
        prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM))
      
      names(ff) <- names(counts)
      
      fit <- MatrixList(lapply(ff, function(f) matrix(f[["pi"]])))
      colnames(fit) <- "null"
      
      lik <- do.call(rbind, lapply(ff, function(f) f[["lik"]])) 
      rownames(lik) <- names(counts)
      colnames(lik) <- "null"
      
      df <- unlist(lapply(ff, function(f) f[["df"]]))
      names(df) <- names(counts)
      
      if(verbose) message("Took ", round(time["elapsed"], 2), " seconds.\n")
      
      return(list(fit = fit, lik = lik, df = df))
      
    })
  
}




bbDS_fitFull_gene <- function(g, counts, 
  ngroups, lgroups, igroups, pi, gamma0, prop_mode, prop_tol, verbose){  
  # g = 1
  if(verbose >= 2)
    message(" Gene:", g)
  
  f <- bb_fitOneGeneManyGroups(y = counts[[g]], ngroups = ngroups, 
    lgroups = lgroups, igroups = igroups, pi = pi[[g]], gamma0 = gamma0[g], 
    verbose = verbose)
  
  return(f)
  
}


bbDS_fitNull_gene <- function(g, counts, 
  pi, gamma0, prop_mode, prop_tol, verbose){  
  # g = 1
  # message("Gene:", g)
  
  f <- bb_fitOneGeneOneGroup(y = counts[[g]], pi = pi[[g]], gamma0 = gamma0[g], 
    verbose = verbose)
  
  return(f)
  
}



bbDS_fitOneModel <- function(counts, samples, pi, dispersion, model = "full",
  verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  inds <-  1:length(counts)
  
  if(length(dispersion) == 1){
    gamma0 <- rep(dispersion, length(inds))
  } else {
    gamma0 <- dispersion
  }
  
  switch(model, 
    
    full = {
      
      if(verbose) message("* Fitting full Beta-Binomial model.. \n")
      
      group <- samples$group
      ngroups <- nlevels(group)
      lgroups <- levels(group)
      
      igroups <- lapply(lgroups, function(gr){which(group == gr)})
      names(igroups) <- lgroups
      
      time <- system.time(ff <- BiocParallel::bplapply(inds, 
        bbDS_fitFull_gene, counts = counts, 
        ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
        pi, gamma0 = gamma0, verbose = verbose, BPPARAM = BPPARAM))
      
      names(ff) <- names(counts)
      
      fit <- MatrixList(lapply(ff, function(f) f[["pi"]]))
      
      lik <- do.call(rbind, lapply(ff, function(f) f[["lik"]])) 
      rownames(lik) <- rownames(counts@unlistData)
      
      df <- rep(NA, nrow(lik))
      names(df) <- rownames(counts@unlistData)
      
      if(verbose >= 2) message("\n")
      if(verbose) message("Took ", round(time["elapsed"], 2), " seconds.\n")
      
      return(list(fit = fit, lik = lik, df = df))
      
    },
    
    null = {
      
      if(verbose) message("* Fitting null Beta-Binomial model.. \n")
      
      time <- system.time(ff <- BiocParallel::bplapply(inds, 
        bbDS_fitNull_gene, counts = counts, 
        pi, gamma0 = gamma0, verbose = verbose, BPPARAM = BPPARAM))

      names(ff) <- names(counts)
      
      fit <- MatrixList(lapply(ff, function(f) matrix(f[["pi"]])))
      colnames(fit) <- "null"
      
      lik <- do.call(rbind, lapply(ff, function(f) matrix(f[["lik"]]))) 
      rownames(lik) <- rownames(counts@unlistData)
      colnames(lik) <- "null"
      
      # df <- unlist(lapply(ff, function(f) f[["df"]]))
      # names(df) <- rownames(counts@unlistData)
      
      df <- rep(1, nrow(lik))
      names(df) <- rownames(counts@unlistData)
      
      if(verbose) message("Took ", round(time["elapsed"], 2), " seconds.\n")
      
      return(list(fit = fit, lik = lik, df = df))
      
    })
  
}


















