dmSQTL_profileLikCommon <- function(disp, counts, genotypes, 
  disp_adjust = TRUE, prop_mode = "constrOptim", prop_tol = 1e-12, 
  verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  if(verbose >= 2) message("Gamma in optimize:", disp)
  
  fit_full <- dmSQTL_fitOneModel(counts = counts, genotypes = genotypes, 
    dispersion = disp, model = "full", prop_mode = prop_mode, 
    prop_tol = prop_tol, verbose=verbose, BPPARAM = BPPARAM)
  
  lik <- sum(unlist(lapply(fit_full, function(g) 
    sum(g@metadata, na.rm = TRUE))), na.rm = TRUE)
  
  if(verbose >= 2) message("lik:", lik)
  
  if(!disp_adjust)
    return(lik)
  
  ## Cox-Reid adjustement
  if(verbose >= 2) message("* Calculating adjustement.. \n")
  
  time <- system.time(adj <- dmSQTL_adjustmentCommon(disp, counts, genotypes, 
    prop = fit_full, BPPARAM = BPPARAM))
  
  if(verbose >= 2) message("Took ", round(time["elapsed"]), " seconds.\n")
  
  adjLik <- lik - adj
  
  if(verbose >= 2) message("adjLik:", adjLik)
  
  return(adjLik)
  
}

#' @importFrom stats optimize

dmSQTL_estimateCommonDispersion <- function(counts, genotypes, 
  disp_adjust = TRUE, 
  disp_interval = c(0, 1e+5), disp_tol = 1e+01, prop_mode = "constrOptim", 
  prop_tol = 1e-12, verbose = FALSE, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  if(verbose) message("* Estimating common dispersion.. \n")
  
  time <- system.time(optimum <- optimize(f = dmSQTL_profileLikCommon, 
    interval = disp_interval, counts = counts, genotypes = genotypes, 
    disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, 
    verbose = max(0, verbose-1), 
    BPPARAM = BPPARAM, maximum = TRUE, tol = disp_tol) )
  
  dispersion <- optimum$maximum
  
  if(verbose) message("Took ", round(time["elapsed"], 2), " seconds.\n")
  
  return(dispersion)
  
}


