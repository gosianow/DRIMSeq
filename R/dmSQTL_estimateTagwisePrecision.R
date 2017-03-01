#' @importFrom stats complete.cases

dmSQTL_estimateTagwisePrecision <- function(counts, genotypes, 
  mean_expression, prec_adjust = TRUE, 
  prec_init = 100, prec_grid_length = 21, prec_grid_range = c(-10, 10), 
  prec_moderation = "none", prec_prior_df = 0, prec_span = 0.1, 
  one_way = TRUE, group_formula = ~ group,
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  time_start <- Sys.time()
  if(verbose) message("* Estimating genewise precision.. \n")
  
  ### Standard grid like in edgeR
  spline_pts <- seq(from = prec_grid_range[1], to = prec_grid_range[2], 
    length = prec_grid_length)
  
  spline_prec <- prec_init * 2^spline_pts
  
  # Calculate the likelihood for each gene and snp 
  # at the spline precision points
  loglik <- matrix(NA, nrow = nrow(genotypes), ncol = prec_grid_length)
  
  for(i in seq(prec_grid_length)){
    # i = 1
    
    loglik[, i] <- dmSQTL_profileLik(prec = spline_prec[i], counts = counts, 
      genotypes = genotypes, prec_adjust = prec_adjust, 
      one_way = one_way, group_formula = group_formula,
      prop_mode = prop_mode, prop_tol = prop_tol, 
      coef_mode = coef_mode, coef_tol = coef_tol,
      verbose = max(0, verbose - 1), BPPARAM = BPPARAM)
    
  }
  
  not_nas <- complete.cases(loglik)        
  loglik <- loglik[not_nas, , drop = FALSE]
  
  if(nrow(loglik) == 0){
    precision <- rep(NA, nrow(genotypes))
    names(precision) <- rownames(genotypes)
    precision <- relist(precision, genotypes@partitioning)
    return(precision)
  }
  
  if(prec_moderation != "none"){
    
    mean_expression <- rep(mean_expression, elementNROWS(genotypes))[not_nas]
    
    loglik <- dm_profileLikModeration(loglik = loglik, 
      mean_expression = mean_expression, 
      prec_moderation = prec_moderation, 
      prec_prior_df = prec_prior_df, prec_span = prec_span)
    
  }
  
  out <- edgeR::maximizeInterpolant(spline_pts, loglik)
  
  # Set NA for genes that tagwise prec could not be calculated 
  precision <- rep(NA, nrow(genotypes))
  names(precision) <- rownames(genotypes)
  precision[not_nas] <- prec_init * 2^out
  precision <- relist(precision, genotypes@partitioning)
  
  time_end <- Sys.time()
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  return(precision)
  
}






