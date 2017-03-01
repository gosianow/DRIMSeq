

dm_CRadjustmentGroups <- function(y, ngroups, lgroups, igroups, 
  gamma0, pi){  
  # pi matrix q x ngroups
  
  if(all(is.na(pi[1, ])))
    return(NA)
  
  adj <- numeric(ngroups)
  
  for(gr in 1:ngroups){
    # gr=1
    
    pi_tmp <- pi[, gr]
    
    if(is.na(pi_tmp[1])){
      
      adj[gr] <- NA
      
    }else{
      
      y_tmp  <- y[, igroups[[gr]], drop = FALSE]
      a <- dm_CRadjustmentOneGroup(y = y_tmp, gamma0, pi = pi_tmp)
      adj[gr] <- a
      
    }
    
  }
  
  adj <- sum(adj, na.rm = TRUE)
  
  if(abs(adj) == Inf)
    return(NA) 
  
  return(adj)
  
}

