# Fitting Dirichlet-multinomial model

dmDS_fitManyGroups_gene <- function(g, counts, 
  ngroups, lgroups, igroups, disp, prop_mode, prop_tol, verbose){  

  if(verbose >= 2)
    message(" Gene:", g)
  
  f <- dm_fitManyGroups(y = counts[[g]], 
    ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
    disp = disp[g], prop_mode = prop_mode, prop_tol = prop_tol, 
    verbose = verbose)
  
  return(f)
  
}

dmDS_fitRegression_gene <- function(g, counts, 
  design, disp, coef_mode, coef_tol, verbose){  

  if(verbose >= 2)
    message(" Gene:", g)
  
  f <- dm_fitRegression(y = counts[[g]], 
    design = design, disp = disp[g], 
    coef_mode = coef_mode, coef_tol = coef_tol, verbose = verbose)
  
  return(f)
  
}

dmDS_fit <- function(counts, design, dispersion,
  one_way = TRUE,
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12, 
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  time_start <- Sys.time()
  if(verbose) message("* Fitting the DM model.. \n")
  
  inds <-  1:length(counts)
  
  # Prepare dispersion
  if(length(dispersion) == 1){
    disp <- rep(dispersion, length(inds))
  } else {
    disp <- dispersion
  }
  
  # Approach from edgeR:
  # If the design is equivalent to a oneway layout, use a shortcut algorithm 
  groups <- edgeR::designAsFactor(design)
  
  if(nlevels(groups) == ncol(design) && one_way){
    
    groups <- factor(groups, labels = paste0("gr", levels(groups)))
    ngroups <- nlevels(groups)
    lgroups <- levels(groups)
    igroups <- lapply(lgroups, function(gr){which(groups == gr)})
    names(igroups) <- lgroups
    
    ff <- BiocParallel::bplapply(inds, dmDS_fitManyGroups_gene, 
      counts = counts, 
      ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
      disp = disp, prop_mode = prop_mode, prop_tol = prop_tol, 
      verbose = verbose, BPPARAM = BPPARAM)
    
    names(ff) <- names(counts)
    
    lik <- unlist(lapply(ff, function(f) sum(f[["lik"]]))) 
    
    prop <- MatrixList(lapply(ff, function(f) f[["prop"]]))
    
    fit <- prop[, groups]
    colnames(fit) <- colnames(counts)
    
    # Get the coefficients like in edgeR::mglmOneWay
    design_unique <- unique(design)
    
    # Use the last feature (qth) as a denominator
    logit_prop <- MatrixList(lapply(ff, function(f) 
      t(t(f[["prop"]])/f[["prop"]][nrow(f[["prop"]]), ])))
    logit_prop <- log(logit_prop@unlistData)
    
    coef <- t(solve(design_unique, t(logit_prop)))
    
    coef <- new("MatrixList", unlistData = coef, 
      partitioning = prop@partitioning)
    
    
  }else{
    
    ff <- BiocParallel::bplapply(inds, dmDS_fitRegression_gene, 
      counts = counts, design = design, disp = disp, 
      coef_mode = coef_mode, coef_tol = coef_tol, 
      verbose = verbose, BPPARAM = BPPARAM)
    
    names(ff) <- names(counts)
    
    lik <- unlist(lapply(ff, function(f) f[["lik"]])) 
    
    coef <- MatrixList(lapply(ff, function(f) f[["b"]]))
    colnames(coef) <- colnames(design)
    
    fit <- MatrixList(lapply(ff, function(f) f[["fit"]]))
    colnames(fit) <- colnames(counts)

  }
  
  time_end <- Sys.time()
  if(verbose >= 2) message("\n")
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  # fit is a MatrixList
  # lik is a vector
  # coef is a MatrixList of matrices q x p
  return(list(fit = fit, lik = lik, coef = coef))
  
  
}


# -----------------------------------------------------------------------------
# Fitting Beta-binomial model
# Currently, recalculating the BB likelihoods using the DM fittings and 
# coefficients

bbDS_fitManyGroups_gene <- function(g, counts, prop,
  ngroups, lgroups, igroups, disp, verbose){  

  if(verbose >= 2)
    message(" Gene:", g)
  
  f <- bb_fitManyGroups(y = counts[[g]], prop = prop[[g]],
    ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
    disp = disp[g], verbose = verbose)
  
  return(f)
  
}



bbDS_fit <- function(counts, fit, design, dispersion,
  one_way = TRUE, 
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  time_start <- Sys.time()
  if(verbose) message("* Fitting the BB model.. \n")
  
  inds <-  1:length(counts)
  
  # Prepare dispersion
  if(length(dispersion) == 1){
    disp <- rep(dispersion, length(inds))
  } else {
    disp <- dispersion
  }
  
  # Approach from edgeR:
  # If the design is equivalent to a oneway layout, use a shortcut algorithm 
  groups <- edgeR::designAsFactor(design)
  
  if(nlevels(groups) == ncol(design) && one_way){
    
    groups <- factor(groups, labels = paste0("gr", levels(groups)))
    ngroups <- nlevels(groups)
    lgroups <- levels(groups)
    igroups <- lapply(lgroups, function(gr){which(groups == gr)})
    names(igroups) <- lgroups
    
    # Use proportions estimated with the DM model
    prop <- fit[, unlist(lapply(igroups, function(x){x[1]}))]
    
    ff <- BiocParallel::bplapply(inds, bbDS_fitManyGroups_gene, 
      counts = counts, prop = prop,
      ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
      disp = disp,
      verbose = verbose, BPPARAM = BPPARAM)
    
    lik <- lapply(ff, function(f){rowSums(f[["lik"]])})
    names(lik) <- NULL
    lik <- unlist(lik)
    names(lik) <- rownames(counts)
    
    # Get the coefficients like in edgeR::mglmOneWay
    design_unique <- unique(design)
    
    logit_prop <- MatrixList(lapply(ff, function(f){
      f[["prop"]]/(1 - f[["prop"]])
    }))
    logit_prop <- log(logit_prop@unlistData)
    
    coef <- t(solve(design_unique, t(logit_prop)))
    
    coef <- new("MatrixList", unlistData = coef, 
      partitioning = prop@partitioning)
    
    
  }else{
    stop("Currently, regression framework is not implemented here!")
  }
  
  time_end <- Sys.time()
  if(verbose >= 2) message("\n")
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  # fit is a MatrixList
  # lik is a vector
  # coef is a MatrixList
  return(list(fit = fit, lik = lik, coef = coef))
  
  
}


















