

dmDS_profileLik <- function(disp, counts, design, 
  disp_adjust = TRUE, prop_mode = "constrOptim", prop_tol = 1e-12, 
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  fit <- dmDS_fit(counts = counts, design = design, dispersion = disp,
  prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, 
  BPPARAM = BPPARAM)
  
  if(!disp_adjust)
    return(fit$lik)
  
  adj <- dmDS_CRadjustment(counts = counts, fit = fit$fit, design = design, 
    dispersion = disp, verbose = verbose, BPPARAM = BPPARAM) 
  
  adj_lik <- fit$lik - adj
  
  return(adj_lik)
  
}

