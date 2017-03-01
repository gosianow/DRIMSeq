

dmSQTL_profileLikCommon <- function(disp, counts, genotypes, 
  disp_adjust = TRUE, disp_interval = c(0, 1e+5), disp_tol = 1e+01,
  one_way = TRUE, group_formula = ~ group,
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  if(verbose >= 2) message("Gamma in optimize:", disp)
  
  adj_lik <- dmSQTL_profileLik(disp = disp, counts = counts, 
    genotypes = genotypes, disp_adjust = disp_adjust, 
    one_way = one_way, group_formula = group_formula,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol,
    verbose = verbose, BPPARAM = BPPARAM)
  
  adj_lik_common <- sum(adj_lik, na.rm = TRUE)
  
  return(adj_lik_common)
  
}

#' @importFrom stats optimize

dmSQTL_estimateCommonDispersion <- function(counts, genotypes, 
  disp_adjust = TRUE, disp_interval = c(0, 1e+5), disp_tol = 1e+01,
  one_way = TRUE, group_formula = ~ group,
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  time_start <- Sys.time()
  if(verbose) message("* Estimating common dispersion.. \n")
  
  optimum <- optimize(f = dmSQTL_profileLikCommon, interval = disp_interval,
    counts = counts, genotypes = genotypes, disp_adjust = disp_adjust,
    one_way = one_way, group_formula = group_formula,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol,
    verbose = max(0, verbose-1), BPPARAM = BPPARAM,
    maximum = TRUE, tol = disp_tol)
  
  dispersion <- optimum$maximum
  
  time_end <- Sys.time()
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  return(dispersion)
  
}


