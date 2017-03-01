

dm_CRadjustmentManyGroups <- function(y, ngroups, lgroups, igroups, 
  disp, prop){  
  # prop matrix q x ngroups
  
  if(all(is.na(prop[1, ])))
    return(NA)
  
  adj <- numeric(ngroups)
  
  for(gr in 1:ngroups){
    # gr=1
    
    prop_tmp <- prop[, gr]
    
    if(is.na(prop_tmp[1])){
      
      adj[gr] <- NA
      
    }else{
      
      y_tmp  <- y[, igroups[[gr]], drop = FALSE]
      a <- dm_CRadjustmentOneGroup(y = y_tmp, disp, prop = prop_tmp)
      adj[gr] <- a
      
    }
    
  }
  
  adj <- sum(adj, na.rm = TRUE)
  
  if(abs(adj) == Inf)
    return(NA) 
  
  return(adj)
  
}

