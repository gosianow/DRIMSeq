##############################################################################
### Tagwise profile likelihood to be optimized
##############################################################################

# gamma0 = splineDisp[i]; y = counts[[g]]

# gamma0 = splineDisp[j]; y = yg

dm_profileLikTagwise <- function(gamma0, y, ngroups, lgroups, igroups, 
  disp_adjust = TRUE, prop_mode = "constrOptimG", prop_tol = 1e-12, 
  verbose = FALSE){
  
  fit <- dm_fitOneGeneManyGroups(y = y, ngroups = ngroups, 
    lgroups = lgroups, igroups = igroups, gamma0 = gamma0, 
    prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose) 
  
  lik <- sum(fit$stats, na.rm = TRUE)
  
  if(!disp_adjust)
    return(lik)
  
  adj <- dm_adjustmentOneGeneManyGroups(y = y, ngroups = ngroups, 
    lgroups = lgroups, igroups = igroups, gamma0 = gamma0, pi = fit$pi) 
  
  adjLik <- lik - adj
  
  return(adjLik)
  
}

