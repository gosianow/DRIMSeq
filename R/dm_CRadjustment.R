
dm_adjustmentOneGeneOneGroup <- function(y, gamma0, pi){
  ### y must be exons vs. samples
  ### If something is wrong, return NAs
  
  ### check for 0s in rows (features)
  keep_row <- rowSums(y) > 0
  y <- y[keep_row, , drop=FALSE]
  
  ### check for 0s in cols (replicates)
  keep_col <- colSums(y) > 0
  y <- y[, keep_col, drop=FALSE]
  
  N <- ncol(y) 
  pi <- pi[keep_row]
  
  adj <- log(det(N * (- dm_HessianG(pi = pi[-length(pi)], gamma0, y)) ))/2 
  ## with Gamma functions ## if pi is NULL then:
  # Error in is.data.frame(x) :
  # dims [product 6] do not match the length of object [0]
  
  return(adj)
  
}


dm_adjustmentOneGeneManyGroups <- function(y, ngroups, lgroups, igroups, 
  gamma0, pi){  
  
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
      a <- dm_adjustmentOneGeneOneGroup(y = y_tmp, gamma0, pi = pi_tmp)
      adj[gr] <- a
      
    }
    
  }
  
  adj <- sum(adj, na.rm = TRUE)
  
  if(abs(adj) == Inf)
    return(NA) 
  
  return(adj)
  
}

