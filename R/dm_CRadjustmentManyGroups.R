

dm_CRadjustmentManyGroups <- function(y, ngroups, lgroups, igroups, 
  prec, prop){  
  # y can not have any rowSums(y) == 0 - assured during dmFilter
  # prop matrix q x ngroups
  
  if(any(is.na(prop[1, ])))
    return(NA)
  
  adj <- numeric(ngroups)
  
  for(gr in 1:ngroups){
    # gr=1
    
    prop_tmp <- prop[, gr]
    y_tmp  <- y[, igroups[[gr]], drop = FALSE]
    
    a <- dm_CRadjustmentOneGroup(y = y_tmp, prec, prop = prop_tmp)
    
    adj[gr] <- a
    
  }
  
  adj <- sum(adj)
  
  if(abs(adj) == Inf)
    return(NA) 
  
  return(adj)
  
}

