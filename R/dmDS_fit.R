

dmDS_fitGroups_gene <- function(g, counts, 
  ngroups, lgroups, igroups, 
  gamma0, prop_mode, prop_tol, verbose){  
  # g = 1
  
  if(verbose >= 2)
    message(" Gene:", g)
  
  f <- dm_fitGroups(y = counts[[g]], 
    ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
    gamma0 = gamma0[g], prop_mode = prop_mode, prop_tol = prop_tol, 
    verbose = verbose)
  
  return(f)
  
}



dmDS_fit <- function(counts, design, dispersion,
  prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, 
  BPPARAM = BiocParallel::SerialParam()){
  
  if(verbose) message("* Fitting the DM model.. \n")
  
  inds <-  1:length(counts)
  
  # Prepare dispersion
  if(length(dispersion) == 1){
    gamma0 <- rep(dispersion, length(inds))
  } else {
    gamma0 <- dispersion
  }
  
  # If the design is equivalent to a oneway layout, use a shortcut algorithm
  groups <- edgeR::designAsFactor(design)
  
  if(nlevels(groups) == ncol(design)){
    
    groups <- factor(groups, labels = paste0("gr", levels(groups)))
    ngroups <- nlevels(groups)
    lgroups <- levels(groups)
    igroups <- lapply(lgroups, function(gr){which(groups == gr)})
    names(igroups) <- lgroups
    
    time <- system.time(ff <- BiocParallel::bplapply(inds, 
      dmDS_fitGroups_gene, counts = counts, 
      ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
      gamma0 = gamma0, prop_mode = prop_mode, prop_tol = prop_tol, 
      verbose = verbose, BPPARAM = BPPARAM))
    
    names(ff) <- names(counts)
    
    lik <- unlist(lapply(ff, function(f) sum(f[["lik"]]))) 
    
    pi <- MatrixList(lapply(ff, function(f) f[["pi"]]))
    
    fit <- pi[, groups]
    colnames(fit) <- colnames(counts)
    
    # Get the coefficients like in edgeR::mglmOneWay
    design_unique <- unique(design)
    
    logit_pi <- MatrixList(lapply(ff, function(f) 
      t(t(f[["pi"]])/f[["pi"]][nrow(f[["pi"]]), ])))
    logit_pi <- log(logit_pi@unlistData)
    
    coeffs <- t(solve(design_unique, t(logit_pi)))
    
    coeffs <- new("MatrixList", unlistData = coeffs, 
      partitioning = pi@partitioning)
    
    
  }else{
    stop("Currently, regression framework is not implemented!")
  }
  
  
  if(verbose >= 2) message("\n")
  if(verbose) message("Took ", round(time["elapsed"], 4), " seconds.\n")
  
  # fit is a MatrixList
  # lik is a vector
  # coeffs is a MatrixList
  return(list(fit = fit, lik = lik, coeffs = coeffs))
  
  
}




bbDS_fitFull_gene <- function(g, counts, 
  ngroups, lgroups, igroups, pi, gamma0, prop_mode, prop_tol, verbose){  
  # g = 1
  if(verbose >= 2)
    message(" Gene:", g)
  
  f <- bb_fitOneGeneManygroupss(y = counts[[g]], ngroups = ngroups, 
    lgroups = lgroups, igroups = igroups, pi = pi[[g]], gamma0 = gamma0[g], 
    verbose = verbose)
  
  return(f)
  
}


bbDS_fitNull_gene <- function(g, counts, 
  pi, gamma0, prop_mode, prop_tol, verbose){  
  # g = 1
  # message("Gene:", g)
  
  f <- bb_fitOneGeneOnegroups(y = counts[[g]], pi = pi[[g]], gamma0 = gamma0[g], 
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
      
      groups <- samples$group
      ngroups <- nlevels(groups)
      lgroups <- levels(groups)
      
      igroups <- lapply(lgroups, function(gr){which(groups == gr)})
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


















