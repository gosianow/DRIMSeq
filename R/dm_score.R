##############################################################################
# Score-function for k-1 parameters (pi) 
##############################################################################

## Using gamma functions
dm_scoreG <- function(pi, gamma0, y){  
  
  k <- nrow(y)
  N <- ncol(y)
  yk <- y[k,]
  y <- y[-k, , drop=FALSE]
  pik <- 1-sum(pi)
  
  S <- gamma0 * rowSums( digamma(y + pi * gamma0) - 
      digamma(pi * gamma0) - matrix(digamma(yk + gamma0 * pik) - 
          digamma(gamma0 * pik), nrow = k-1, ncol = N, byrow = TRUE) ) 
  
  return(S)
  
} 



dm_score <- function(pi, gamma0, y){  
  
  k <- nrow(y)
  N <- ncol(y) 
  yk <- y[k, ]
  y <- y[-k, , drop=FALSE]
  pik <- 1-sum(pi)
  
  S <- rep(0, k-1)
  
  for(j in 1:N){ 
    # j=1
    if(yk[j] == 0) Skj <- 0
    else Skj <- sum(gamma0 / (pik * gamma0 + 1:yk[j] - 1))    
    
    for(i in 1:(k-1)){
      # i=1
      if(y[i,j] == 0) Sij <- 0
      else Sij <- sum(gamma0 / (pi[i] * gamma0 + 1:y[i,j] - 1)) 
      S[i] <- S[i] + Sij
    }
    S <- S - Skj
  }
  
  return(S)
  
}





dm_score_regG <- function(b, x, gamma0, y){  
  
  y <- t(y)
  
  n <- nrow(x)
  q <- ncol(y)
  p <- ncol(x)
  
  b <- matrix(b, p, q-1)
  
  z <- exp(x %*% b)
  pi <- z/(1 + rowSums(z))
  piq <- 1 - rowSums(pi)
  
  yq <- y[, q]
  y <- y[, -q, drop = FALSE]
  
  
  S <- t(x) %*% ((digamma(y + pi*gamma0) - 
      digamma(pi*gamma0) + digamma(yq + piq*gamma0) - 
      digamma(piq*gamma0)) * pi * (1 - pi) * gamma0)
  
  
  return(as.numeric(S))
  
  
  
} 








