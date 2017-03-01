##############################################################################
## Computes the log-likelihood -- with gamma functions -- for q-1 parameters
##############################################################################

# pi <- pi[-length(pi)]

dm_likG <- function(pi, gamma0, y){
  # pi has length of q-1
  # gamma0 has length 1
  # y has q rows and n number of columns
  # This function returns likelihhod without normalizing component, 
  # but it is OK for the LR test
  
  n <- ncol(y)
  m_i <- colSums(y)
  
  pi <- c(pi, 1 - sum(pi))
  
  l <- n * lgamma(gamma0) - sum(lgamma(m_i + gamma0)) + 
    sum( colSums( lgamma(y + pi * gamma0) - lgamma(pi * gamma0) ) )
  
  # normalizing_part <- sum(lgamma(m_i + 1) - colSums(lgamma(y + 1)))
  
  return(l)
  
}


bb_likG <- function(pi, gamma0, y){
  # pi has length of q
  # gamma0 has length 1
  # y has q rows and n number of columns
 
 m_i <- colSums(y)
 
 l <- rep(NA, length(pi))
 
 for(i in 1:length(pi)){
   
  l[i] <- dm_likG(pi[i], gamma0, y = rbind(y[i, ], m_i - y[i, ]))
  
 }
  
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
# log-likelihood for q-1 parameters (pi)
##############################################################################

# pi <- pi[-length(pi)]

dm_lik <- function(pi, gamma0, y){
  # pi has length q-1
  # gamma0 has length 1
  # y has q rows and n number of columns
  
  q <- nrow(y)
  n <- ncol(y)
  m_i <- colSums(y)  
  l <- 0
  
  pi <- c(pi, 1 - sum(pi))
  
  for(i in 1:n){  
    # i=1
    l <- l - sum(log(gamma0 + 1:m_i[i] - 1))   
    for(j in 1:q){   
      # j=3
      if(y[j,i] == 0) lji <- 0
      else lji <- sum(log(pi[j] * gamma0 + 1:y[j,i] - 1))     
      l <- l + lji      
    }
  }
  
  return(l)
  
}

