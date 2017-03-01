##############################################################################
## estimate pi for given dispersion
##############################################################################

#' @importFrom stats constrOptim

dm_fitOneGeneOneGroup <- function(y, gamma0, prop_mode = c("constrOptim", 
  "constrOptimG")[2], prop_tol = 1e-12, verbose = FALSE){
  ### y must be features vs. samples
  ### If something is wrong, return NAs
  
  # NAs for genes with one feature
  kk <- nrow(y)
  if(kk < 2 || is.na(gamma0)) 
    return(list(pi = rep(NA, kk), lik = NA, df = NA))
  
  ### check for 0s in rows (features)
  keep_row <- rowSums(y) > 0
  ### must be at least two features
  if(sum(keep_row) < 2) 
    return(list(pi = rep(NA, kk), lik = NA, df = NA))
  
  y <- y[keep_row, , drop=FALSE]
  
  ### check for 0s in columns (replicates)
  keep_col <- colSums(y) > 0
  y <- y[, keep_col, drop=FALSE]
  
  pi_init <- rowSums(y)/sum(y)
  k <- length(pi_init) ## k - number of features
  
  if(sum(keep_col) == 1){
    
    keep_row[keep_row] <- pi_init
    pi <- keep_row
    
    df <- k - 1
    lik <- dm_likG(pi = pi_init[-k], gamma0 = gamma0, y = y)
    
    return(list(pi = pi, lik = lik, df = df))
  }
  
  switch(prop_mode, 
    
    ### must have constraint for SUM pi = 1 --> 
    ### sum(pi) < 1 + eps & sum(pi) > 1 - eps
    constrOptim = { ## for k-1 parameters
      # if(verbose) message("\n gene:", colnames(y)[1], "gamma0:", gamma0)
      
      ui <- rbind(diag(rep(1, k-1), k-1), diag(rep(-1, k-1), k-1), rep(-1, k-1))
      ci <- c(rep(0, k-1), rep(-1, k-1), -1 + .Machine$double.eps) 
      # ui <- rbind(diag(rep(1, k-1)), diag(rep(-1, k-1)))
      # ci <- c(rep(0, k-1), rep(-1, k-1))
      
      co <- constrOptim(pi_init[-k], f = dm_lik, grad = dm_score, 
        ui = ui, ci = ci, control = list(fnscale = -1, reltol = prop_tol), 
        gamma0 = gamma0, y = y)
      
      pi <- co$par
      pi <- c(pi, 1-sum(pi))
      lik <- co$value
      
    }, 
    
    constrOptimG = { ## for k-1 parameters with Gamma functions
      # if(verbose) message("\n gene:", colnames(y)[1], "gamma0:", gamma0)
      
      ui <- rbind(diag(rep(1, k-1), k-1), diag(rep(-1, k-1), k-1), rep(-1, k-1))
      ci <- c(rep(0, k-1), rep(-1, k-1), -1 + .Machine$double.eps) 
      # ui <- rbind(diag(rep(1, k-1)), diag(rep(-1, k-1)))
      # ci <- c(rep(0, k-1), rep(-1, k-1))
      
      co <- constrOptim(pi_init[-k], f = dm_likG, grad = dm_scoreG, 
        ui = ui, ci = ci, control = list(fnscale = -1, reltol = prop_tol), 
        gamma0 = gamma0, y = y)
      
      pi <- co$par
      pi <- c(pi, 1-sum(pi))
      lik <- co$value
      
    })
  
  keep_row[keep_row] <- pi
  pi <- keep_row
  
  df <- k - 1
  
  return(list(pi = pi, lik = lik, df = df))
  
}



bb_fitOneGeneOneGroup <- function(y, pi, gamma0, verbose = FALSE){
  
  # NAs for genes with one feature
  kk <- nrow(y)
  if(kk < 2 || is.na(gamma0)) 
    return(list(pi = rep(NA, kk), lik = NA, df = NA))
  
  ### check for 0s in rows (features)
  keep_row <- rowSums(y) > 0
  ### must be at least two features
  if(sum(keep_row) < 2) 
    return(list(pi = rep(NA, kk), lik = NA, df = NA))
  
  y <- y[keep_row, , drop=FALSE]
  pi <- pi[keep_row]
  
  ### check for 0s in columns (replicates)
  keep_col <- colSums(y) > 0
  y <- y[, keep_col, drop=FALSE]
  
  lik <- rep(NA, length(keep_row))
  
  lik[keep_row] <- bb_likG(pi = pi, gamma0 = gamma0, y = y)
    
  keep_row[keep_row] <- pi
  pi <- keep_row
  
  df <- rep(1, length(keep_row))
  
  return(list(pi = pi, lik = lik, df = df))
  
}




