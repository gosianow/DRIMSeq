##############################################################################
## Computes the deviance -- with gamma functions -- for q-1 parameters
##############################################################################

dm_devG <- function(prop, prec, y){
  ## prop has length of q-1
  ## prec has length 1
  ## y has q rows and n columns 
  
  prop <- c(prop, 1 - sum(prop))
  
  ll_mod <- sum(lgamma(y + prop * prec) - lgamma(prop * prec) , na.rm = TRUE )
  
  prop_sat <- y/matrix(colSums(y), nrow(y), ncol(y), byrow = TRUE)
  
  ll_sat <- sum(lgamma(y + prop_sat * prec) - lgamma(prop_sat * prec) , 
    na.rm = TRUE) # Inf for y = 0
  
  D <- 2 * (ll_sat - ll_mod)
  
  return(D)
  
}
