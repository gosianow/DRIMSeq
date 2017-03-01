##############################################################################
## estimate pi for given dispersion
##############################################################################

#' @importFrom stats constrOptim

dm_fitOneGroup <- function(y, gamma0, prop_mode = c("constrOptim", 
  "constrOptimG")[2], prop_tol = 1e-12, verbose = FALSE){
  ### y must be features vs. samples
  ### If something is wrong, return NAs
  
  q <- nrow(y)
  
  # NAs for genes with one feature
  if(q < 2 || is.na(gamma0)) 
    return(list(pi = rep(NA, q), lik = NA))
  
  # Check for 0s in rows (features)
  keep_row <- rowSums(y) > 0
  # Must be at least two features
  if(sum(keep_row) < 2) 
    return(list(pi = rep(NA, q), lik = NA))
    
  # Last feature can not be zero since
  # we use the last feature as a denominator in logit
  if(keep_row[q] == 0)
    return(list(pi = rep(NA, q), lik = NA))

  y <- y[keep_row, , drop=FALSE]
  q <- nrow(y)
  
  # Check for 0s in columns (replicates)
  keep_col <- colSums(y) > 0
  y <- y[, keep_col, drop=FALSE]
  
  pi_init <- rowSums(y)/sum(y)
  
  if(sum(keep_col) == 1){
    
    lik <- dm_likG(pi = pi_init[-q], gamma0 = gamma0, y = y)
    
    keep_row[keep_row] <- pi_init
    pi <- keep_row
    
    return(list(pi = pi, lik = lik))
  }
  
  switch(prop_mode, 
    
    ### Must have constraint for SUM pi = 1 --> 
    ### sum(pi) < 1 + eps & sum(pi) > 1 - eps
    constrOptim = { 
    ## For q-1 parameters
      
      ui <- rbind(diag(rep(1, q-1), q-1), diag(rep(-1, q-1), q-1), rep(-1, q-1))
      ci <- c(rep(0, q-1), rep(-1, q-1), -1 + .Machine$double.eps) 
      # ui <- rbind(diag(rep(1, q-1)), diag(rep(-1, q-1)))
      # ci <- c(rep(0, q-1), rep(-1, q-1))
      
      co <- constrOptim(pi_init[-q], f = dm_lik, grad = dm_score, 
        ui = ui, ci = ci, control = list(fnscale = -1, reltol = prop_tol), 
        gamma0 = gamma0, y = y)
      
      pi <- co$par
      pi <- c(pi, 1-sum(pi))
      lik <- co$value
      
    }, 
    
    constrOptimG = { 
    # For q-1 parameters with Gamma functions
      
      ui <- rbind(diag(rep(1, q-1), q-1), diag(rep(-1, q-1), q-1), rep(-1, q-1))
      ci <- c(rep(0, q-1), rep(-1, q-1), -1 + .Machine$double.eps) 
      # ui <- rbind(diag(rep(1, q-1)), diag(rep(-1, q-1)))
      # ci <- c(rep(0, q-1), rep(-1, q-1))
      
      co <- constrOptim(pi_init[-q], f = dm_likG, grad = dm_scoreG, 
        ui = ui, ci = ci, control = list(fnscale = -1, reltol = prop_tol), 
        gamma0 = gamma0, y = y)
      
      pi <- co$par
      pi <- c(pi, 1-sum(pi))
      lik <- co$value
      
    })
  
  keep_row[keep_row] <- pi
  pi <- keep_row
  
  # pi numeric vector of length q
  # lik numeric of lenght 1
  return(list(pi = pi, lik = lik))
  
}



bb_fitOneGroup <- function(y, pi, gamma0, verbose = FALSE){
  # Recalculates likelihood for BB, where pi is estimated with DM
  
  q <- nrow(y)
  
  # BB lik only for non-NA pi from DM
  if(any(is.na(pi)))
    return(list(pi = rep(NA, q), lik = rep(NA, q)))
  
  # NAs for genes with one feature
  if(q < 2 || is.na(gamma0)) 
    return(list(pi = rep(NA, q), lik = rep(NA, q)))
  
  # Check for 0s in rows (features with zero proportions)
  keep_row <- rowSums(y) > 0
  # Must be at least two non zero features
  if(sum(keep_row) < 2) 
    return(list(pi = rep(NA, q), lik = rep(NA, q)))
  
  y <- y[keep_row, , drop=FALSE]
  pi <- pi[keep_row]
  
  # Check for 0s in columns (replicates)
  keep_col <- colSums(y) > 0
  y <- y[, keep_col, drop=FALSE]
  
  lik <- rep(NA, q)
  
  lik[keep_row] <- bb_likG(pi = pi, gamma0 = gamma0, y = y)
  
  keep_row[keep_row] <- pi
  pi <- keep_row

  # pi numeric vector of length q
  # lik numeric vector of length q
  return(list(pi = pi, lik = lik))
  
}




















