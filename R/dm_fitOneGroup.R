##############################################################################
## estimate prop for given dispersion
##############################################################################

#' @importFrom stats constrOptim

dm_fitOneGroup <- function(y, disp, 
  prop_mode = "constrOptim", prop_tol = 1e-12){
  # y matrix q x n
  # If something is wrong, return NAs
  
  q <- nrow(y)
  
  # NAs for genes with one feature
  if(q < 2 || is.na(disp)) 
    return(list(prop = rep(NA, q), lik = NA))
  
  # Check for 0s in rows (features)
  keep_row <- rowSums(y) > 0
  # There must be at least two non-zero features
  if(sum(keep_row) < 2) 
    return(list(prop = rep(NA, q), lik = NA))
  
  # Last feature can not be zero since
  # we use the last feature as a denominator in logit
  if(keep_row[q] == 0)
    return(list(prop = rep(NA, q), lik = NA))
  
  y <- y[keep_row, , drop=FALSE]
  q <- nrow(y)
  
  # Check for 0s in columns (replicates)
  keep_col <- colSums(y) > 0
  y <- y[, keep_col, drop=FALSE]
  
  prop_init <- rowSums(y)/sum(y)
  
  # If there is only one replicate, use empirical props as output
  if(sum(keep_col) == 1){
    
    lik <- dm_likG(prop = prop_init[-q], disp = disp, y = y)
    
    keep_row[keep_row] <- prop_init
    prop <- keep_row
    
    return(list(prop = prop, lik = lik))
  }
  
  switch(prop_mode, 
    
    constrOptim = { 
      ### Must have constraint for SUM prop = 1 --> 
      ### sum(prop) < 1 + eps & sum(prop) > 1 - eps 
      
      ui <- rbind(diag(rep(1, q-1), q-1), diag(rep(-1, q-1), q-1), rep(-1, q-1))
      ci <- c(rep(0, q-1), rep(-1, q-1), -1 + .Machine$double.eps) 
      # ui <- rbind(diag(rep(1, q-1)), diag(rep(-1, q-1)))
      # ci <- c(rep(0, q-1), rep(-1, q-1))
      
      # Minimization
      co <- constrOptim(theta = prop_init[-q], f = dm_likG_neg, 
        grad = dm_scoreG_neg, 
        ui = ui, ci = ci, control = list(reltol = prop_tol), 
        method = "BFGS",
        disp = disp, y = y)
      
      if(co$convergence == 0){
        prop <- co$par
        prop <- c(prop, 1-sum(prop))
        lik <- -co$value
      }else{
        return(list(prop = rep(NA, length(keep_row)), lik = NA))
      }
      
    })
  
  keep_row[keep_row] <- prop
  prop <- keep_row
  
  # prop numeric vector of length q
  # lik numeric of lenght 1
  return(list(prop = prop, lik = lik))
  
}



bb_fitOneGroup <- function(y, disp, prop){
  # Recalculates likelihood for BB, where prop is estimated with DM
  
  q <- nrow(y)
  
  # BB lik only for non-NA prop from DM
  if(any(is.na(prop)))
    return(list(prop = rep(NA, q), lik = rep(NA, q)))
  
  # NAs for genes with one feature
  if(q < 2 || is.na(disp)) 
    return(list(prop = rep(NA, q), lik = rep(NA, q)))
  
  # Check for 0s in rows (features with zero proportions)
  keep_row <- rowSums(y) > 0
  # Must be at least two non zero features
  if(sum(keep_row) < 2) 
    return(list(prop = rep(NA, q), lik = rep(NA, q)))
  
  y <- y[keep_row, , drop=FALSE]
  prop <- prop[keep_row]
  
  # Check for 0s in columns (replicates)
  keep_col <- colSums(y) > 0
  y <- y[, keep_col, drop=FALSE]
  
  lik <- rep(NA, q)
  
  lik[keep_row] <- bb_likG(prop = prop, disp = disp, y = y)
  
  keep_row[keep_row] <- prop
  prop <- keep_row
  
  # prop numeric vector of length q
  # lik numeric vector of length q
  return(list(prop = prop, lik = lik))
  
}




















