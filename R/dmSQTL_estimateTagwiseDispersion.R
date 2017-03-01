##############################################################################
# calculate tagwise dispersions 
##############################################################################

#' @importFrom stats complete.cases


dmSQTL_grid_dm_profileLikTagwise <- function(g, counts, genotypes, 
  disp_grid_length, seq_disp_grid_length, splineDisp, disp_adjust, 
  prop_mode, prop_tol, verbose){
  # g = 1
  
  y <- counts[[g]]
  snps <- genotypes[[g]]
  ll <- matrix(0, nrow(snps), disp_grid_length)
  
  
  for(i in 1:nrow(snps)){
    # i = 1
    
    NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
    yg <- y[, !NAs]             
    group <- snps[i, !NAs]
    group <- factor(group)
    ngroups <- nlevels(group)
    lgroups <- levels(group)
    nlibs <- length(group)
    
    igroups <- lapply(lgroups, function(gr){which(group == gr)})
    names(igroups) <- lgroups
    
    for(j in seq_disp_grid_length){
      # j = 1 
      
      out <- dm_profileLikTagwise(disp = splineDisp[j], y = yg, 
        ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
        disp_adjust = disp_adjust, prop_mode = prop_mode, 
        prop_tol = prop_tol, verbose = verbose)
      
      if(is.na(out)){
        ll[i, ] <- NA
        break
      }
      
      ll[i, j] <- out
      
    }
  }
  
  return(ll)
  
}



dmSQTL_estimateTagwiseDispersion <- function(counts, genotypes, 
  mean_expression, disp_adjust = TRUE, disp_mode = "grid",
  disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = 100, 
  disp_init_weirMoM = TRUE, disp_grid_length = 21, 
  disp_grid_range = c(-10, 10),  disp_moderation = "none", disp_prior_df = 0, 
  disp_span = 0.1, prop_mode = "constrOptim", prop_tol = 1e-12, 
  verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  inds <- 1:length(counts)
  if(verbose) message("* Estimating genewise dispersion.. \n")
  
  
        ### Standard grid like in edgeR
        splinePts <- seq(from = disp_grid_range[1], to = disp_grid_range[2], 
          length = disp_grid_length)
        splineDisp <- disp_init * 2^splinePts
        seq_disp_grid_length <- seq(disp_grid_length)
        
        ### Calculate the likelihood for each gene at the spline dispersion points
        loglikL <- BiocParallel::bplapply(inds, 
          dmSQTL_grid_dm_profileLikTagwise, 
          counts = counts, genotypes = genotypes, 
          disp_grid_length = disp_grid_length, 
          seq_disp_grid_length = seq_disp_grid_length, 
          splineDisp = splineDisp, disp_adjust = disp_adjust, 
          prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, 
          BPPARAM = BPPARAM)
        
        loglik <- do.call(rbind, loglikL)
        not_nas <- complete.cases(loglik)        
        loglik <- loglik[not_nas, , drop = FALSE]
        
        if(disp_moderation != "none"){
          
          mean_expression <- rep(mean_expression, width(genotypes))[not_nas]
          
          loglik <- dm_profileLikModeration(loglik = loglik, mean_expression = mean_expression, disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, disp_span = disp_span)
          
        }
        
        out <- edgeR::maximizeInterpolant(splinePts, loglik)
        
        #### set NA for genes that tagwise disp could not be calculated            
        dispersion <- rep(NA, length(not_nas))
        names(dispersion) <- rownames(genotypes@unlistData)
        dispersion[not_nas] <- disp_init * 2^out
        dispersion <- relist(dispersion, genotypes@partitioning)
        
  
  if(verbose) message("Took ", round(time["elapsed"], 2), " seconds.\n")
  
  return(dispersion)
}






