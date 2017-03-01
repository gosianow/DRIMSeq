
dm_fitGroups <- function(y, ngroups, lgroups, igroups, 
  gamma0, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE) {
  
  q <- nrow(y)
  
  pi <- matrix(NA, nrow = q, ncol = ngroups, dimnames = list(rownames(y), 
    lgroups))
  lik <- rep(NA, ngroups)
  names(lik) <- lgroups
  
  if(is.na(gamma0) || q < 2) 
    return(list(pi = pi, lik = lik))
  
  
  for(gr in 1:ngroups){
    # gr = 1
    
    fit_gr <- dm_fitOneGroup(y = y[, igroups[[gr]], 
      drop = FALSE], gamma0 = gamma0, prop_mode = prop_mode, 
      prop_tol = prop_tol, verbose = verbose)
    
    
    if(is.na(fit_gr[["lik"]])) {
      pi <- matrix(NA, nrow = q, ncol = ngroups, dimnames = list(rownames(y), 
        lgroups))
      lik <- rep(NA, ngroups)
      names(lik) <- lgroups
      return(list(pi = pi, lik = lik))
    }
    
    
    pi[, gr] <- fit_gr[["pi"]] 
    lik[gr] <- fit_gr[["lik"]]
    
  }
  
  # pi and lik can have NAs
  # pi matrix q x ngroups
  # lik vector of length ngroups
  return(list(pi = pi, lik = lik))  
  
}


bb_fitGroups <- function(y, ngroups, lgroups, igroups, 
  pi, gamma0, verbose = FALSE){
  # This function calculates BB likelihoods 
  # Proportions pi are estimated with DM model
  
  q <- nrow(y)
  
  lik <- matrix(NA, nrow = q, ncol = ngroups, dimnames = list(rownames(y), 
    lgroups))

  
  if(is.na(gamma0) || q < 2) 
    return(list(pi = pi, lik = lik))
  
  
  for(gr in 1:ngroups){
    # gr = 1
    
    fit_gr <- bb_fitOneGroup(y = y[, igroups[[gr]], drop = FALSE], 
      gamma0 = gamma0, pi = pi[, gr], verbose = verbose)
    
    lik[, gr] <- fit_gr[["lik"]]
    
  }
  
  lik[rowSums(is.na(lik)) > 0, ] <- rep(NA, ngroups)
  
  # pi and lik can have NAs
  # pi matrix q x ngroups
  # lik vector of length ngroups
  return(list(pi = pi, lik = lik)) 
  
}















