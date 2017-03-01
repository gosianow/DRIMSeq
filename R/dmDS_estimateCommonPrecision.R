#' @importFrom stats optimize


dmDS_profileLikCommon <- function(prec, counts, design, 
  prec_adjust = TRUE, one_way = TRUE, 
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  if(verbose >= 2) message("Gamma in optimize:", prec)
  
  adj_lik <- dmDS_profileLik(prec = prec, counts = counts, design = design, 
    prec_adjust = prec_adjust, one_way = one_way, 
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol,
    verbose = verbose, BPPARAM = BPPARAM)
  
  adj_lik_common <- sum(adj_lik, na.rm = TRUE)
  
  return(adj_lik_common)
  
}


dmDS_estimateCommonPrecision <- function(counts, design, 
  prec_adjust = TRUE, prec_interval = c(0, 1e+5), prec_tol = 1e-01, 
  one_way = TRUE,
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  time_start <- Sys.time()
  if(verbose) message("* Estimating common precision.. \n")
  
  optimum <- optimize(f = dmDS_profileLikCommon, 
    interval = prec_interval,
    counts = counts, design = design, 
    prec_adjust = prec_adjust, one_way = one_way,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol,
    verbose = max(0, verbose-1), BPPARAM = BPPARAM,
    maximum = TRUE, tol = prec_tol)
  
  precision <- optimum$maximum
  
  time_end <- Sys.time()
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  return(precision)
  
}

