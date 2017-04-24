# Fitting the Dirichlet-multinomial model

dmSQTL_fitManyGroups_gene <- function(g, counts, genotypes,
  prec, prop_mode, prop_tol, verbose){  
  # g = 6
  
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
      prec = prec[[g]][i], prop_mode = prop_mode, prop_tol = prop_tol)
    
    lik <- sum(f$lik)
    
    fit <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
    fit[, !NAs] <- f$prop[, groups]
    
    return(list(lik = lik, fit = fit))
    
  })
  
  lik <- unlist(lapply(ff, function(f) f[["lik"]])) 
  fit <- MatrixList(lapply(ff, function(f) f[["fit"]]))
  colnames(fit) <- colnames(counts)
  
  # lik is a vector of length nrow(x)
  # fit is a MatrixList of matrices q x n
  # DOES NOT return coef
  return(list(lik = lik, fit = fit))
  
}

#' @importFrom stats model.matrix

dmSQTL_fitRegression_gene <- function(g, counts, genotypes, 
  group_formula = ~ group,
  prec, coef_mode, coef_tol, verbose){  
  
  if(verbose >= 2) message(" Gene:", g)
  
  y <- counts[[g]]
  x <- genotypes[[g]]
  
  ff <- lapply(1:nrow(x), function(i){
    # i = 2
    
    NAs <- is.na(x[i, ]) | is.na(y[1, ])            
    yy <- y[, !NAs, drop = FALSE]             
    xx <- x[i, !NAs]
    
    design <- model.matrix(group_formula, data = data.frame(group = xx))
    
    f <- dm_fitRegression(y = yy, 
      design = design, prec = prec[[g]][i], 
      coef_mode = coef_mode, coef_tol = coef_tol)
    
    fit <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
    fit[, !NAs] <- f$fit
    
    return(list(lik = f$lik, fit = fit))
    
  })
  
  lik <- unlist(lapply(ff, function(f) f[["lik"]])) 
  fit <- MatrixList(lapply(ff, function(f) f[["fit"]]))
  colnames(fit) <- colnames(counts)
  
  # lik is a vector of length nrow(x)
  # fit is a MatrixList of matrices q x n
  # DOES NOT return coef
  return(list(lik = lik, fit = fit))
  
}


dmSQTL_fit <- function(counts, genotypes, precision,
  one_way = TRUE, group_formula = ~ group,
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12, 
  return_fit = FALSE, return_coef = FALSE, 
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  time_start <- Sys.time()
  if(verbose) message("* Fitting the DM model.. \n")
  
  inds <-  1:length(counts)
  
  # Prepare precision
  if(class(precision) == "numeric"){ 
    prec <- relist(rep(precision, nrow(genotypes)), genotypes@partitioning)
  } else {
    prec <- precision
  }
  
  # Approach from edgeR glmFit.default:
  # If oneway layout, use a shortcut algorithm 
  if(one_way){
    
    if(verbose) message("   Using the one way approach. \n")
    
    ff <- BiocParallel::bplapply(inds, dmSQTL_fitManyGroups_gene, 
      counts = counts, genotypes = genotypes,
      prec = prec, prop_mode = prop_mode, prop_tol = prop_tol, 
      verbose = verbose, BPPARAM = BPPARAM)
    
    names(ff) <- names(counts)
    
    lik <- lapply(ff, function(f) f[["lik"]])
    
    if(return_fit){
      fit <- lapply(ff, function(f) f[["fit"]])
    }else{
      fit <- list()  
    }
    
  }else{
    
    if(verbose) message("   Using the regression approach. \n")
    
    ff <- BiocParallel::bplapply(inds, dmSQTL_fitRegression_gene, 
      counts = counts, genotypes = genotypes, group_formula = group_formula,
      prec = prec, coef_mode = coef_mode, coef_tol = coef_tol, 
      verbose = verbose, BPPARAM = BPPARAM)
    
    names(ff) <- names(counts)
    
    lik <- lapply(ff, function(f) f[["lik"]])
    
    if(return_fit){
      fit <- lapply(ff, function(f) f[["fit"]])
    }else{
      fit <- list()  
    }
    
  }
  
  time_end <- Sys.time()
  if(verbose >= 2) message("\n")
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  # fit is a list of length G of MatrixLists
  # lik is a list of length G ofvectors
  # coef Currently, do not compute coef 
  return(list(fit = fit, lik = lik, coef = list()))
  
}


















