##############################################################################
## Computes the log-likelihood -- with gamma functions -- for k-1 parameters
##############################################################################

# pi <- pi[-length(pi)]

dm_likG <- function(pi, gamma0, y){
  ## pi has length of k-1
  ## gamma0 has length 1
  ## y has k rows and any number of columns n
  ## This function returns likelihhod without normalizing component, 
  ## but it is OK for LR test.
  
  N <- ncol(y)
  S <- colSums(y)
  
  pi <- c(pi, 1 - sum(pi))
  
  l <- sum( colSums( lgamma(y + pi * gamma0) - lgamma(pi * gamma0) ) )
  
  l <- N * lgamma(gamma0) - sum(lgamma(S + gamma0)) + l
  
  # normalizing_part <- sum(lgamma(S + 1) - colSums(lgamma(y + 1)))
  
  return(l)
  
}


dm_lik_regG <- function(b, x, gamma0, y){
  ## b has length of p * (q-1)
  ## gamma0 has length 1
  ## y has q rows and n columns
  ## This function returns likelihhod without normalizing component, 
  ## but it is OK for LR test.
  
  y <- t(y)
  
  n <- nrow(x)
  q <- ncol(y)
  p <- ncol(x)
  
  b <- matrix(b, p, q-1)
  
  z <- exp(x %*% b)
  pi <- z/(1 + rowSums(z))
  pi <- cbind(pi, 1 - rowSums(pi))
  
  s <- rowSums(y)
  
  l <- sum(rowSums(lgamma(y + gamma0 * pi) - lgamma(pi * gamma0)))
  
  l <- n * lgamma(gamma0) - sum(lgamma(s + gamma0)) + l
  
  # normalizing_part <- sum(lgamma(s + 1) - rowSums(lgamma(y + 1)))
  
  return(l)
  
}

##############################################################################
# log-likelihood for k-1 parameters (pi)
##############################################################################

# pi <- pi[-length(pi)]

dm_lik <- function(pi, gamma0, y){
  ## pi has length of k-1
  ## gamma0 has length 1
  ## y has k rows and any number of columns
  
  k <- nrow(y)
  N <- ncol(y)
  S <- colSums(y)  
  l <- 0
  
  pi <- c(pi, 1 - sum(pi))
  
  for(j in 1:N){  
    # j=1
    l <- l - sum(log(gamma0 + 1:S[j] - 1))   
    for(i in 1:k){   
      # i=3
      if(y[i,j] == 0) lij <- 0
      else lij <- sum(log(pi[i] * gamma0 + 1:y[i,j] - 1))     
      l <- l + lij      
    }
  }
  
  return(l)
  
}

