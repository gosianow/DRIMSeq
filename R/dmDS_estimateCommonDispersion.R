#' @importFrom stats optimize


dmDS_profileLikCommon <- function(gamma0, counts, samples, disp_adjust = TRUE, 
  prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  if(verbose >= 2) message("Gamma in optimize:", gamma0)
  
  fit <- dmDS_fitOneModel(counts = counts, samples = samples, 
    dispersion = gamma0, model = "full", prop_mode = prop_mode, 
    prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  
  lik <- sum(fit[["lik"]], na.rm = TRUE) 
  
  if(verbose >= 2) message("lik:", lik)
  
  if(!disp_adjust)
    return(lik)
  
  ## Cox-Reid adjustement for common dispersion
  if(verbose >= 2) message("* Calculating adjustement.. \n")
  
  time <- system.time(adj <- dmDS_adjustmentCommon(gamma0, counts = counts, 
    samples = samples, pi = fit[["fit"]], BPPARAM = BPPARAM))
  
  if(verbose >= 2) message("Took ", round(time["elapsed"]), " seconds.\n")
  
  adjLik <- lik - adj
  
  if(verbose >= 2) message("adjLik:", adjLik)
  
  
  return(adjLik)
  
}


dmDS_estimateCommonDispersion <- function(counts, samples, disp_adjust = TRUE, 
  disp_interval = c(0, 1e+5), disp_tol = 1e-01, prop_mode = "constrOptimG", 
  prop_tol = 1e-12, verbose = FALSE, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  if(verbose) message("* Estimating common dispersion.. \n")
  
  time <- system.time(optimum <- optimize(f = dmDS_profileLikCommon, 
    interval = disp_interval,
    counts = counts, samples = samples, disp_adjust = disp_adjust, 
    prop_mode = prop_mode, prop_tol = prop_tol, 
    verbose = max(0, verbose-1), BPPARAM = BPPARAM,
    maximum = TRUE, tol = disp_tol) )
  
  dispersion <- optimum$maximum
  
  if(verbose) message("Took ", round(time["elapsed"], 2), " seconds.\n")
  
  return(dispersion)
  
}

