

dmDS_profileLik <- function(prec, counts, design, 
  prec_adjust = TRUE, one_way = TRUE, 
  prop_mode = "constrOptim", prop_tol = 1e-12,
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  fit <- dmDS_fit(counts = counts, design = design, precision = prec,
    one_way = one_way,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol,
    verbose = verbose, BPPARAM = BPPARAM)
  
  if(!prec_adjust)
    return(fit$lik)
  
  adj <- dmDS_CRadjustment(counts = counts, fit = fit$fit, design = design, 
    precision = prec, one_way = one_way,
    verbose = verbose, BPPARAM = BPPARAM) 
  
  adj_lik <- fit$lik - adj
  
  # adj_lik has length G
  return(adj_lik)
  
}

