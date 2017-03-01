#' @importFrom stats optimize


dmDS_profileLikCommon <- function(disp, counts, design, disp_adjust = TRUE, 
  prop_mode = "constrOptim", prop_tol = 1e-12, verbose = FALSE, 
  BPPARAM = BiocParallel::SerialParam()){
  
  if(verbose >= 2) message("Gamma in optimize:", disp)
  
  adj_lik <- dmDS_profileLik(disp = disp, counts = counts, design = design, 
    disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, 
    verbose = verbose, BPPARAM = BPPARAM)
  
  adj_lik_common <- sum(adj_lik, na.rm = TRUE)
  
  return(adj_lik_common)
  
}


dmDS_estimateCommonDispersion <- function(counts, design, disp_adjust = TRUE, 
  disp_interval = c(0, 1e+5), disp_tol = 1e-01, prop_mode = "constrOptim", 
  prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  time_start <- Sys.time()
  if(verbose) message("* Estimating common dispersion.. \n")
  
  optimum <- optimize(f = dmDS_profileLikCommon, 
    interval = disp_interval,
    counts = counts, design = design, disp_adjust = disp_adjust, 
    prop_mode = prop_mode, prop_tol = prop_tol, 
    verbose = max(0, verbose-1), BPPARAM = BPPARAM,
    maximum = TRUE, tol = disp_tol)
  
  dispersion <- optimum$maximum
  
  time_end <- Sys.time()
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  return(dispersion)
  
}

