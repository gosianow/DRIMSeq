

dmSQTL_profileLikCommon <- function(prec, counts, genotypes, 
  prec_adjust = TRUE, prec_interval = c(0, 1e+5), prec_tol = 1e+01,
  one_way = TRUE, group_formula = ~ group,
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  if(verbose >= 2) message("Gamma in optimize:", prec, "\n")
  
  adj_lik <- dmSQTL_profileLik(prec = prec, counts = counts, 
    genotypes = genotypes, prec_adjust = prec_adjust, 
    one_way = one_way, group_formula = group_formula,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol,
    verbose = verbose, BPPARAM = BPPARAM)
  
  adj_lik_common <- sum(adj_lik, na.rm = TRUE)
  
  return(adj_lik_common)
  
}

#' @importFrom stats optimize

dmSQTL_estimateCommonPrecision <- function(counts, genotypes, 
  prec_adjust = TRUE, prec_interval = c(0, 1e+5), prec_tol = 1e+01,
  one_way = TRUE, group_formula = ~ group,
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  time_start <- Sys.time()
  if(verbose) message("* Estimating common precision.. \n")
  
  optimum <- optimize(f = dmSQTL_profileLikCommon, interval = prec_interval,
    counts = counts, genotypes = genotypes, prec_adjust = prec_adjust,
    one_way = one_way, group_formula = group_formula,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol,
    verbose = max(0, verbose-1), BPPARAM = BPPARAM,
    maximum = TRUE, tol = prec_tol)
  
  precision <- optimum$maximum
  
  time_end <- Sys.time()
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  return(precision)
  
}


