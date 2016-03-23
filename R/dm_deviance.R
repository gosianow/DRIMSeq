##############################################################################
## Computes the deviance -- with gamma functions -- for k-1 proportion parameters
##############################################################################

# pi <- pi[-length(pi)]

dm_devG <- function(pi, gamma0, y){
  ## pi has length of k-1
  ## gamma0 has length 1
  ## y has k rows and any number of columns n
  
  pi <- c(pi, 1 - sum(pi))
  
  ll_mod <- sum(lgamma(y + pi * gamma0) - lgamma(pi * gamma0) , na.rm = TRUE )
  
  pi_sat <- y/matrix(colSums(y), nrow(y), ncol(y), byrow = TRUE)
  
  ll_sat <- sum(lgamma(y + pi_sat * gamma0) - lgamma(pi_sat * gamma0) , 
    na.rm = TRUE) # Inf for y = 0
  
  D <- 2 * (ll_sat - ll_mod)
  
  return(D)
  
}
