# Fitting Dirichlet-multinomial model

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
  
  time_start <- Sys.time()
  if(verbose) message("* Fitting the DM model.. \n")
  
  inds <-  1:length(counts)
  
  # Prepare dispersion
  if(length(dispersion) == 1){
    gamma0 <- rep(dispersion, length(inds))
  } else {
    gamma0 <- dispersion
  }
  
  # Approach from edgeR:
  # If the design is equivalent to a oneway layout, use a shortcut algorithm 
  groups <- edgeR::designAsFactor(design)
  
  if(nlevels(groups) == ncol(design)){
    
    groups <- factor(groups, labels = paste0("gr", levels(groups)))
    ngroups <- nlevels(groups)
    lgroups <- levels(groups)
    igroups <- lapply(lgroups, function(gr){which(groups == gr)})
    names(igroups) <- lgroups
    
    ff <- BiocParallel::bplapply(inds, dmDS_fitGroups_gene, 
      counts = counts, 
      ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
      gamma0 = gamma0, prop_mode = prop_mode, prop_tol = prop_tol, 
      verbose = verbose, BPPARAM = BPPARAM)
    
    names(ff) <- names(counts)
    
    lik <- unlist(lapply(ff, function(f) sum(f[["lik"]]))) 
    
    pi <- MatrixList(lapply(ff, function(f) f[["pi"]]))
    
    fit <- pi[, groups]
    colnames(fit) <- colnames(counts)
    
    # Get the coefficients like in edgeR::mglmOneWay
    design_unique <- unique(design)
    
    # Use the last feature (qth) as a denominator
    logit_pi <- MatrixList(lapply(ff, function(f) 
      t(t(f[["pi"]])/f[["pi"]][nrow(f[["pi"]]), ])))
    logit_pi <- log(logit_pi@unlistData)
    
    coeffs <- t(solve(design_unique, t(logit_pi)))
    
    coeffs <- new("MatrixList", unlistData = coeffs, 
      partitioning = pi@partitioning)
    
    
  }else{
    stop("Currently, regression framework is not implemented!")
  }
  
  time_end <- Sys.time()
  if(verbose >= 2) message("\n")
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  # fit is a MatrixList
  # lik is a vector
  # coeffs is a MatrixList
  return(list(fit = fit, lik = lik, coeffs = coeffs))
  
  
}


# -----------------------------------------------------------------------------
# Fitting Beta-binomial model
# Currently, recalculating the BB likelihood using the DM fitting and 
# coefficients

bbDS_fitGroups_gene <- function(g, counts, pi,
  ngroups, lgroups, igroups, 
  gamma0, verbose){  
  # g = 1
  
  if(verbose >= 2)
    message(" Gene:", g)
  
  f <- bb_fitGroups(y = counts[[g]], pi = pi[[g]],
    ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
    gamma0 = gamma0[g], verbose = verbose)
  
  return(f)
  
}



bbDS_fit <- function(counts, fit, design, dispersion,
  prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, 
  BPPARAM = BiocParallel::SerialParam()){
  
  time_start <- Sys.time()
  if(verbose) message("* Fitting the BB model.. \n")
  
  inds <-  1:length(counts)
  
  # Prepare dispersion
  if(length(dispersion) == 1){
    gamma0 <- rep(dispersion, length(inds))
  } else {
    gamma0 <- dispersion
  }
  
  # Approach from edgeR:
  # If the design is equivalent to a oneway layout, use a shortcut algorithm 
  groups <- edgeR::designAsFactor(design)
  
  if(nlevels(groups) == ncol(design)){
    
    groups <- factor(groups, labels = paste0("gr", levels(groups)))
    ngroups <- nlevels(groups)
    lgroups <- levels(groups)
    igroups <- lapply(lgroups, function(gr){which(groups == gr)})
    names(igroups) <- lgroups
    
    # Use proportions estimated with the DM model
    pi <- fit[, unlist(lapply(igroups, function(x){x[1]}))]
    
    ff <- BiocParallel::bplapply(inds, bbDS_fitGroups_gene, 
      counts = counts, pi = pi,
      ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
      gamma0 = gamma0,
      verbose = verbose, BPPARAM = BPPARAM)
    
    names(ff) <- names(counts)
    
    lik <- lapply(ff, function(f){rowSums(f[["lik"]])})
    names(lik) <- NULL
    lik <- unlist(lik)
    
    # Get the coefficients like in edgeR::mglmOneWay
    design_unique <- unique(design)
    
    logit_pi <- MatrixList(lapply(ff, function(f){
      f[["pi"]]/(1 - f[["pi"]])
      }))
    logit_pi <- log(logit_pi@unlistData)
    
    coeffs <- t(solve(design_unique, t(logit_pi)))
    
    coeffs <- new("MatrixList", unlistData = coeffs, 
      partitioning = pi@partitioning)
    
    
  }else{
    stop("Currently, regression framework is not implemented!")
  }
  
  time_end <- Sys.time()
  if(verbose >= 2) message("\n")
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  # fit is a MatrixList
  # lik is a vector
  # coeffs is a MatrixList
  return(list(fit = fit, lik = lik, coeffs = coeffs))
  
  
}


















