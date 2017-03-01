
dm_fitOneGeneManyGroups <- function(y, ngroups, lgroups, igroups, 
  gamma0, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE) {
  
  k <- nrow(y)
  
  pi <- matrix(NA, nrow = k, ncol = ngroups, dimnames = list(rownames(y), 
    lgroups))
  lik <- rep(NA, ngroups)
  names(lik) <- lgroups
  
  if (is.na(gamma0) || k < 2) 
    return(list(pi = pi, lik = lik))
  
  
  for (gr in 1:ngroups) {
    # gr = 1
    
    fit_gr <- dm_fitOneGeneOneGroup(y = y[, igroups[[gr]], 
      drop = FALSE], gamma0 = gamma0, prop_mode = prop_mode, 
      prop_tol = prop_tol, verbose = verbose)
    
    
    if (is.na(fit_gr[["lik"]])) {
      pi <- matrix(NA, nrow = k, ncol = ngroups, dimnames = list(rownames(y), 
        lgroups))
      lik <- rep(NA, ngroups)
      names(lik) <- lgroups
      return(list(pi = pi, lik = lik))
    }
    
    
    pi[, gr] <- fit_gr[["pi"]] 
    lik[gr] <- fit_gr[["lik"]]
    
  }
  
  return(list(pi = pi, lik = lik))  ### pi and lik can have NAs
  
}


