

dmSQTL_profileLik <- function(disp, counts, genotypes, 
  disp_adjust = TRUE, one_way = TRUE, group_formula = ~ group,
  prop_mode = "constrOptim", prop_tol = 1e-12,
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  fit <- dmSQTL_fit(counts = counts, genotypes = genotypes, dispersion = disp,
    one_way = one_way, group_formula = group_formula,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol,
    return_fit = TRUE, return_coef = FALSE,
    verbose = verbose, BPPARAM = BPPARAM)
  
  if(!disp_adjust)
    return(unlist(fit$lik))
  
  adj <- dmSQTL_CRadjustment(counts = counts, fit = fit$fit, 
    genotypes = genotypes, group_formula = group_formula,
    dispersion = disp, one_way = one_way,
    verbose = verbose, BPPARAM = BPPARAM) 
  
  adj_lik <- unlist(fit$lik) - unlist(adj)
  
  # adj_lik is a vector containing all the genes and snps
  return(adj_lik)
  
}

