
dm_fitManyGroups <- function(y, ngroups, lgroups, igroups, 
  prec, prop_mode = "constrOptim", prop_tol = 1e-12){
  # y can not have any rowSums(y) == 0 - assured during dmFilter
  
  q <- nrow(y)
  
  prop <- matrix(NA, nrow = q, ncol = ngroups)
  colnames(prop) <- lgroups
  rownames(prop) <- rownames(y)
  lik <- rep(NA, ngroups)
  names(lik) <- lgroups
  
  if(q < 2 || is.na(prec)) 
    return(list(prop = prop, lik = lik))
  
  for(gr in 1:ngroups){
    # gr = 1
    
    fit_gr <- dm_fitOneGroup(y = y[, igroups[[gr]], 
      drop = FALSE], prec = prec, prop_mode = prop_mode, 
      prop_tol = prop_tol) 
    
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


bb_fitManyGroups <- function(y, ngroups, lgroups, igroups, 
  prec, prop){
  # This function calculates BB likelihoods 
  # Proportions prop are estimated with the DM model
  # y can not have any rowSums(y) == 0 - assured during dmFilter
  
  q <- nrow(y)
  
  lik <- matrix(NA, nrow = q, ncol = ngroups)
  colnames(lik) <- lgroups
  
  if(is.na(prec)) 
    return(list(prop = prop, lik = lik))
  
  for(gr in 1:ngroups){
    # gr = 1
    
    fit_gr <- bb_fitOneGroup(y = y[, igroups[[gr]], drop = FALSE], 
      prec = prec, prop = prop[, gr])
    
    lik[, gr] <- fit_gr[["lik"]]
    
  }
  
  lik[rowSums(is.na(lik)) > 0, ] <- NA
  
  # prop and lik can have NAs
  # prop matrix q x ngroups
  # lik matrix q x ngroups
  return(list(prop = prop, lik = lik)) 
  
}















