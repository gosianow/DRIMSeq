##############################################################################
# calculate tagwise dispersions 
##############################################################################

#' @importFrom stats optimize optim constrOptim complete.cases

dmDS_optimize_dm_profileLikTagwise <- function(g, disp_interval, counts, 
  ngroups, lgroups, igroups,
  disp_adjust, prop_mode, prop_tol, verbose, disp_tol){
  # g = 1
  
  ### return NA if gene has 1 exon or observations in one sample in group 
  ### (anyway this gene would not be fitted by dmFit)
  gamma0 <- disp_interval[1] + 
    (1 - (sqrt(5) - 1)/2) * (disp_interval[2] - disp_interval[1])
  
  if(is.na(dm_profileLikTagwise(gamma0 = gamma0, y = counts[[g]], 
    ngroups = ngroups, lgroups = lgroups, igroups = igroups,
    disp_adjust = disp_adjust, prop_mode = prop_mode, 
    prop_tol = prop_tol, verbose = verbose)))
    return(NA) 
  
  optimum <- optimize(f = dm_profileLikTagwise, 
    interval = disp_interval, 
    y = counts[[g]], 
    ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
    disp_adjust = disp_adjust, prop_mode = prop_mode, 
    prop_tol = prop_tol, verbose = verbose,  
    maximum = TRUE, tol = disp_tol) 
  
  return(optimum$maximum)  
  
}


dmDS_optim_dm_profileLikTagwise <- function(g, disp_interval, counts,  
  ngroups, lgroups, igroups,
  disp_init, disp_init_weirMoM, disp_adjust, prop_mode, 
  prop_tol, verbose, disp_tol){
  # g = 12
  
  ### return NA if gene has 1 exon or observations in one sample in group 
  ### (anyway this gene would not be fitted by dmFit)
  gamma0 <- disp_interval[1] + 
    (1 - (sqrt(5) - 1)/2) * (disp_interval[2] - disp_interval[1])
  
  if(is.na(dm_profileLikTagwise(gamma0 = gamma0, y = counts[[g]], 
    ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
    disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, 
    verbose = verbose)))
    return(NA) 
  
  if(disp_init_weirMoM){
    disp_init_tmp <- dm_weirMoM(y = counts[[g]], se = FALSE)
    if(is.na(disp_init_tmp))
      disp_init_tmp <- disp_init
  }else{
    disp_init_tmp <- disp_init
  }
  
  optimum <- NA
  
  try(optimum <- optim(par = disp_init_tmp, fn = dm_profileLikTagwise, 
    gr = NULL, y = counts[[g]], ngroups=ngroups, lgroups=lgroups, 
    igroups=igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, 
    prop_tol = prop_tol, verbose = verbose,  method = "L-BFGS-B", 
    lower = 1e-2, upper = 1e+10, 
    control = list(fnscale = -1, factr = disp_tol))$par , silent = TRUE)
  
  return(optimum)  
  
}


dmDS_constrOptim_dm_profileLikTagwise <- function(g, disp_interval, counts,  
  ngroups, lgroups, igroups,
  disp_init, disp_init_weirMoM, disp_adjust, prop_mode, 
  prop_tol, verbose, disp_tol){
  # g = 1
  
  ### return NA if gene has 1 exon or observations in one sample in group 
  ### (anyway this gene would not be fitted by dmFit)
  gamma0 <- disp_interval[1] + 
    (1 - (sqrt(5) - 1)/2) * (disp_interval[2] - disp_interval[1])
  
  if(is.na(dm_profileLikTagwise(gamma0 = gamma0, y = counts[[g]], 
    ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
    disp_adjust = disp_adjust, prop_mode = prop_mode, 
    prop_tol = prop_tol, verbose = verbose)))
    return(NA) 
  
  ui <- 1
  ci <- 1e-8
  
  if(disp_init_weirMoM){
    disp_init_tmp <- dm_weirMoM(y = counts[[g]], se = FALSE)
    if(is.na(disp_init_tmp))
      disp_init_tmp <- disp_init
  }else{
    disp_init_tmp <- disp_init
  }
  
  optimum <- constrOptim(theta = disp_init_tmp, dm_profileLikTagwise, 
    grad = NULL, method = "Nelder-Mead", ui = ui, ci = ci, 
    control = list(fnscale = -1, reltol = disp_tol), 
    y = counts[[g]], ngroups = ngroups, lgroups = lgroups, 
    igroups = igroups, disp_adjust = disp_adjust, 
    prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
  
  return(optimum$par) 
  
}


dmDS_grid_dm_profileLikTagwise <- function(g, seq_disp_grid_length, splineDisp, 
  counts,
  ngroups, lgroups, igroups,
  disp_adjust, prop_mode, prop_tol, verbose,
  disp_grid_length){
  # g = 1237
  
  ll <- numeric(disp_grid_length)
  
  for(i in seq_disp_grid_length){
    # i = 1
    
    out <- dm_profileLikTagwise(gamma0 = splineDisp[i], y = counts[[g]], 
      ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
      disp_adjust = disp_adjust, prop_mode = prop_mode, 
      prop_tol = prop_tol, verbose = verbose)
    
    if(is.na(out)){
      ll <- rep(NA, disp_grid_length)
      break
    }
    
    ll[i] <- out
    
  }
  
  return(ll)
  
}



dmDS_estimateTagwiseDispersion <- function(counts, samples, mean_expression, 
  disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+5), 
  disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, 
  disp_grid_length = 21, disp_grid_range = c(-10, 10), 
  disp_moderation = "none", disp_prior_df = 10, disp_span = 0.3, 
  prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  inds <- 1:length(counts)
  
  group <- samples$group
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  nlibs <- length(group)
  
  igroups <- lapply(lgroups, function(gr){which(group == gr)})
  names(igroups) <- lgroups
  
  if(verbose) message("* Estimating genewise dispersion.. \n")
  
  time <- system.time(
    
    switch(disp_mode, 
      
      optimize={
        
        disp_list <- BiocParallel::bplapply(inds, 
          dmDS_optimize_dm_profileLikTagwise,
          disp_interval = disp_interval, counts = counts,
          ngroups = ngroups, lgroups = lgroups, igroups = igroups,
          disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, 
          verbose = verbose, disp_tol = disp_tol, 
          BPPARAM = BPPARAM)
        
        names(disp_list) <- names(counts)  
        dispersion <- unlist(disp_list)
        
      },
      
      optim={
        
        disp_list <- BiocParallel::bplapply(inds, 
          dmDS_optim_dm_profileLikTagwise, 
          disp_interval = disp_interval, counts = counts,  
          ngroups = ngroups, lgroups = lgroups, igroups = igroups,
          disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, 
          disp_adjust = disp_adjust, 
          prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, 
          disp_tol = disp_tol, 
          BPPARAM = BPPARAM)
        
        names(disp_list) <- names(counts)  
        dispersion <- unlist(disp_list)
        
      },
      
      constrOptim={
        
        disp_list <- BiocParallel::bplapply(inds, 
          dmDS_constrOptim_dm_profileLikTagwise, 
          disp_interval = disp_interval, counts = counts,  
          ngroups = ngroups, lgroups = lgroups, igroups = igroups,
          disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, 
          disp_adjust = disp_adjust, 
          prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose, 
          disp_tol = disp_tol, 
          BPPARAM = BPPARAM)
        
        names(disp_list) <- names(counts)  
        dispersion <- unlist(disp_list)
        
      },
      
      grid={
        
        splinePts <- seq(from = disp_grid_range[1], to = disp_grid_range[2], 
          length = disp_grid_length)
        splineDisp <- disp_init * 2^splinePts
        
        ### calculate the likelihood for each gene at the spline dispersion points
        seq_disp_grid_length <- seq(disp_grid_length)
        
        loglikL <- BiocParallel::bplapply(inds, dmDS_grid_dm_profileLikTagwise, 
          seq_disp_grid_length = seq_disp_grid_length, splineDisp = splineDisp,
          counts = counts, 
          ngroups = ngroups, lgroups = lgroups, igroups = igroups,
          disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, 
          verbose = verbose, disp_grid_length = disp_grid_length, 
          BPPARAM = BPPARAM)
        
        loglik <- do.call(rbind, loglikL)
        NAs <- complete.cases(loglik)        
        loglik <- loglik[NAs, , drop = FALSE]
        
        if(nrow(loglik) == 0){
          dispersion <- rep(NA, length(inds))
          names(dispersion) <- names(counts)
          return(dispersion)
        }
        
        if(disp_moderation != "none"){
          
          nlibs <- length(group)
          
          ### analogy to edgeR
          # priorN <- disp_prior_df/(nlibs - ngroups) 
          priorN <- disp_prior_df
          
          switch(disp_moderation, 
            
            common={
              
              moderation <- colMeans(loglik)
              loglik <- sweep(loglik, 2, priorN * moderation, FUN = "+")
              
            },
            
            trended={
              
              o <- order(mean_expression[NAs])
              oo <- order(o)
              width <- floor(disp_span * nrow(loglik))
              
              moderation <- edgeR::movingAverageByCol(loglik[o,], 
                width = width)[oo,]
              
              ### like in edataR estimateTagwiseDisp
              loglik <- loglik + priorN * moderation 
              ### like in edataR dispCoxReidInterpolateTagwise
              # loglik <- (loglik + priorN * moderation)/(1 + priorN) 
              
            }
          )
        }
        
        out <- edgeR::maximizeInterpolant(splinePts, loglik)
        
        ### set NA for genes that tagwise disp could not be calculated 
        dispersion <- rep(NA, length(inds))
        names(dispersion) <- names(counts)
        dispersion[NAs] <- disp_init * 2^out
        
      }))
  
  if(verbose) message("Took ", round(time["elapsed"], 2), " seconds.\n")
  
  return(dispersion)
  
}






