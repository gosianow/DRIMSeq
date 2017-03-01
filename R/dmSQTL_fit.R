# Fitting the Dirichlet-multinomial model

dmSQTL_fitManyGroups_gene <- function(g, counts, genotypes,
  disp, prop_mode, prop_tol, verbose){  
  
  if(verbose >= 2) message(" Gene:", g)
  
  y <- counts[[g]]
  x <- genotypes[[g]]
  
  ff <- lapply(1:nrow(x), function(i){
    # i = 2
    
    NAs <- is.na(x[i, ]) | is.na(y[1, ])            
    yy <- y[, !NAs, drop = FALSE]             
    xx <- x[i, !NAs]
    
    groups <- factor(xx)
    ngroups <- nlevels(groups)
    lgroups <- levels(groups)
    igroups <- lapply(lgroups, function(gr){which(groups == gr)})
    names(igroups) <- lgroups
    
    f <- dm_fitManyGroups(y = yy, 
      ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
      disp = disp[[g]][i], prop_mode = prop_mode, prop_tol = prop_tol)
    
    lik <- sum(f$lik)
    
    fit <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
    fit[, !NAs] <- f$prop[, groups]
    
    return(list(lik = lik, fit = fit))
    
  })
  
  lik <- unlist(lapply(ff, function(f) f[["lik"]])) 
  fit <- MatrixList(lapply(ff, function(f) f[["fit"]]))
  
  # lik is a vector of length nrow(x)
  # fit is a MatrixList of matrices q x n
  # DOES NOT return coef
  return(list(lik = lik, fit = fit))
  
}

#' @importFrom stats model.matrix

dmSQTL_fitRegression_gene <- function(g, counts, genotypes, 
  disp, coef_mode, coef_tol, verbose){  
  
  if(verbose >= 2) message(" Gene:", g)
  
  y <- counts[[g]]
  x <- genotypes[[g]]
  
  ff <- lapply(1:nrow(x), function(i){
    # i = 1
    
    NAs <- is.na(x[i, ]) | is.na(y[1, ])            
    yy <- y[, !NAs, drop = FALSE]             
    xx <- x[i, !NAs]
    
    design <- model.matrix(~ group, data = data.frame(group = xx))
    
    f <- dm_fitRegression(y = yy, 
      design = design, disp = disp[[g]][i], 
      coef_mode = coef_mode, coef_tol = coef_tol)
    
    fit <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
    fit[, !NAs] <- f$fit
    
    return(list(lik = f$lik, fit = fit))
    
  })
  
  lik <- unlist(lapply(ff, function(f) f[["lik"]])) 
  fit <- MatrixList(lapply(ff, function(f) f[["fit"]]))
  
  # lik is a vector of length nrow(x)
  # fit is a MatrixList of matrices q x n
  # DOES NOT return coef
  return(list(lik = lik, fit = fit))
  
}





dmSQTL_fit <- function(counts, design, dispersion,
  one_way = TRUE,
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12, 
  return_fit = FALSE, return_coef = FALSE, 
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
    
    if(verbose) message("   Using the one way approach. \n")
    
    groups <- factor(groups, labels = paste0("gr", levels(groups)))
    ngroups <- nlevels(groups)
    lgroups <- levels(groups)
    igroups <- lapply(lgroups, function(gr){which(groups == gr)})
    names(igroups) <- lgroups
    
    ff <- BiocParallel::bplapply(inds, dmSQTL_fitManyGroups_gene, 
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
    
    # Use the last feature (q-th) as a denominator
    logit_prop <- MatrixList(lapply(ff, function(f) 
      t(t(f[["prop"]])/f[["prop"]][nrow(f[["prop"]]), ])))
    logit_prop <- log(logit_prop@unlistData)
    
    coef <- t(solve(design_unique, t(logit_prop)))
    
    coef <- new("MatrixList", unlistData = coef, 
      partitioning = prop@partitioning)
    colnames(coef) <- colnames(design)
    
  }else{
    
    if(verbose) message("   Using the regression approach. \n")
    
    ff <- BiocParallel::bplapply(inds, dmSQTL_fitRegression_gene, 
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
  
  # fit is a MatrixList of matrices q x n
  # lik is a vector of length G
  # coef is a MatrixList of matrices q x p
  return(list(fit = fit, lik = lik, coef = coef))
  
  
}


# -----------------------------------------------------------------------------
# Fitting the Beta-binomial model
# Currently, recalculating the BB likelihoods and coefficients using the 
# DM fittings/proportions

bbSQTL_fitManyGroups_gene <- function(g, counts, prop,
  ngroups, lgroups, igroups, disp, verbose){  
  
  if(verbose >= 2)
    message(" Gene:", g)
  
  f <- bb_fitManyGroups(y = counts[[g]], prop = prop[[g]],
    ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
    disp = disp[g])
  
  return(f)
  
}


bbSQTL_fitRegression_gene <- function(g, counts, 
  design, disp, fit, verbose){  
  
  if(verbose >= 2)
    message(" Gene:", g)
  
  f <- bb_fitRegression(y = counts[[g]], 
    design = design, disp = disp[g], fit = fit[[g]])
  
  return(f)
  
}


bbSQTL_fit <- function(counts, fit, design, dispersion,
  one_way = TRUE, verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
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
    
    if(verbose) message("   Using the one way approach. \n")
    
    groups <- factor(groups, labels = paste0("gr", levels(groups)))
    ngroups <- nlevels(groups)
    lgroups <- levels(groups)
    igroups <- lapply(lgroups, function(gr){which(groups == gr)})
    names(igroups) <- lgroups
    
    # Use proportions estimated with the DM model
    prop <- fit[, unlist(lapply(igroups, function(x){x[1]}))]
    
    # Recalculate BB likelihoods
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
    
    # design_unique must be squared for solve()
    coef <- t(solve(design_unique, t(logit_prop)))
    
    coef <- new("MatrixList", unlistData = coef, 
      partitioning = prop@partitioning)
    
    
  }else{
    
    if(verbose) message("   Using the regression approach. \n")
    
    ff <- BiocParallel::bplapply(inds, bbDS_fitRegression_gene, 
      counts = counts, design = design, disp = disp, 
      fit = fit, verbose = verbose, BPPARAM = BPPARAM)
    
    names(ff) <- names(counts)
    
    lik <- unlist(lapply(ff, function(f) f[["lik"]])) 
    names(lik) <- rownames(counts)
    
    coef <- MatrixList(lapply(ff, function(f) f[["b"]]))
    colnames(coef) <- colnames(design)
    
    fit <- MatrixList(lapply(ff, function(f) f[["fit"]]))
    colnames(fit) <- colnames(counts)
    
  }
  
  time_end <- Sys.time()
  if(verbose >= 2) message("\n")
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  # fit is a MatrixList of matrices q x p
  # lik is a vector of length nrow(counts) = total number of features
  # coef is a MatrixList of matrices q x p
  return(list(fit = fit, lik = lik, coef = coef))
  
  
}


















