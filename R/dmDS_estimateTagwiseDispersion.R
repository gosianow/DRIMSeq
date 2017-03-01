#' @importFrom stats complete.cases


dmDS_estimateTagwiseDispersion <- function(counts, design, mean_expression, 
  disp_adjust = TRUE, 
  disp_init = 100, disp_grid_length = 21, disp_grid_range = c(-10, 10), 
  disp_moderation = "none", disp_prior_df = 0, disp_span = 0.1, 
  one_way = TRUE,
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  time_start <- Sys.time()
  if(verbose) message("* Estimating genewise dispersion.. \n")
  
  ### Standard grid like in edgeR
  spline_pts <- seq(from = disp_grid_range[1], to = disp_grid_range[2], 
    length = disp_grid_length)
  
  spline_disp <- disp_init * 2^spline_pts
  
  # Calculate the likelihood for each gene at the spline dispersion points
  loglik <- matrix(NA, nrow = length(counts), ncol = disp_grid_length)
  
  for(i in seq(disp_grid_length)){
    # i = 1
    
    loglik[, i] <- dmDS_profileLik(disp = spline_disp[i], counts = counts, 
      design = design, disp_adjust = disp_adjust, one_way = one_way,
      prop_mode = prop_mode, prop_tol = prop_tol, 
      coef_mode = coef_mode, coef_tol = coef_tol,
      verbose = max(0, verbose - 1), BPPARAM = BPPARAM)

  }
  
  not_nas <- complete.cases(loglik)        
  loglik <- loglik[not_nas, , drop = FALSE]
  
  if(nrow(loglik) == 0){
    dispersion <- rep(NA, length(counts))
    names(dispersion) <- names(counts)
    return(dispersion)
  }
  
  if(disp_moderation != "none"){
    
    mean_expression <- mean_expression[not_nas]
    
    loglik <- dm_profileLikModeration(loglik = loglik, 
      mean_expression = mean_expression, 
      disp_moderation = disp_moderation, 
      disp_prior_df = disp_prior_df, disp_span = disp_span)
    
  }
  
  out <- edgeR::maximizeInterpolant(spline_pts, loglik)
  
  # Set NA for genes that tagwise disp could not be calculated 
  dispersion <- rep(NA, length(counts))
  names(dispersion) <- names(counts)
  dispersion[not_nas] <- disp_init * 2^out
  
  time_end <- Sys.time()
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  return(dispersion)
  
}






