
dm_fitManyGroups <- function(y, ngroups, lgroups, igroups, 
  disp, prop_mode = "constrOptim", prop_tol = 1e-12, verbose = FALSE) {
  # y can not have any rowSums(y) == 0
  
  q <- nrow(y)
  
  prop <- matrix(NA, nrow = q, ncol = ngroups)
  colnames(prop) <- lgroups
  rownames(prop) <- rownames(y)
  lik <- rep(NA, ngroups)
  names(lik) <- lgroups
  
  if(q < 2 || is.na(disp)) 
    return(list(prop = prop, lik = lik))
  
  for(gr in 1:ngroups){
    # gr = 1
    
    fit_gr <- dm_fitOneGroup(y = y[, igroups[[gr]], 
      drop = FALSE], disp = disp, prop_mode = prop_mode, 
      prop_tol = prop_tol, verbose = verbose) 
    
    if(is.na(fit_gr[["lik"]])) {
      prop <- matrix(NA, nrow = q, ncol = ngroups)
      colnames(prop) <- lgroups
      rownames(prop) <- rownames(y)
      lik <- rep(NA, ngroups)
      names(lik) <- lgroups
      return(list(prop = prop, lik = lik))
    }
    
    prop[, gr] <- fit_gr[["prop"]] 
    lik[gr] <- fit_gr[["lik"]]
    
  }
  
  # prop and lik can have NAs
  # prop matrix q x ngroups
  # lik vector of length ngroups
  return(list(prop = prop, lik = lik))  
  
}


bb_fitManyGroups <- function(y, prop, ngroups, lgroups, igroups, 
  disp, verbose = FALSE){
  # This function calculates BB likelihoods 
  # Proportions prop are estimated with the DM model
  # y can not have any rowSums(y) == 0
  
  q <- nrow(y)
  
  lik <- rep(NA, ngroups)
  names(lik) <- lgroups
  
  if(is.na(disp)) 
    return(list(prop = prop, lik = lik))
  
  for(gr in 1:ngroups){
    # gr = 1
    
    fit_gr <- bb_fitOneGroup(y = y[, igroups[[gr]], drop = FALSE], 
      disp = disp, prop = prop[, gr], verbose = verbose)
    
    lik[, gr] <- fit_gr[["lik"]]
    
  }
  
  lik[rowSums(is.na(lik)) > 0, ] <- rep(NA, ngroups)
  
  # prop and lik can have NAs
  # prop matrix q x ngroups
  # lik vector of length ngroups
  return(list(prop = prop, lik = lik)) 
  
}















