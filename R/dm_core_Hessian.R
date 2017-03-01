##############################################################################
# Hessian for Dirichlet-multinomial model gamma+ and k-1 proportions 
##############################################################################

# pi <- pi[-length(pi)]

### Using gamma functions
dm_HessianG <- function(pi, gamma0, y){  
  ## pi has length of k-1
  
  k <- nrow(y)
  N <- ncol(y)
  D <- rep(0, k-1)
  Dil <- 0
  pik <- 1-sum(pi)
  ykm1 <- y[-k, , drop=FALSE]
  
  Dil <- gamma0^2 * sum(trigamma(y[k,] + gamma0 * pik) - trigamma(gamma0 * pik))
  
  D <- gamma0^2 * rowSums(trigamma(ykm1 + pi * gamma0) - trigamma(pi * gamma0))
  
  H <- matrix(Dil, k-1, k-1)
  diag(H) <- diag(H) + D
  
  return(H)
  
  
}


dm_Hessian <- function(pi, gamma0, y){  
  ## pi has length of k-1
  
  k <- nrow(y)
  N <- ncol(y)
  D <- rep(0, k-1)
  Dil <- 0
  pik <- 1-sum(pi)
  
  for(j in 1:N){
    # j=1
    if(y[k, j] == 0) Dil <- Dil + 0
    else Dil <- Dil +  sum(-gamma0^2 / (pik * gamma0 + 1:y[k, j] - 1) ^2) 
    
    for(i in 1:(k-1)){
      # i=1
      if(y[i,j] == 0) Dii <- 0
      else Dii <- sum(-gamma0^2 / (pi[i] * gamma0 + 1:y[i,j] - 1) ^2) 
      D[i] <- D[i] + Dii
    }
  }
  
  H <- matrix(Dil, k-1, k-1)
  diag(H) <- diag(H) + D
  
  return(H)
  
}




