##############################################################################
## Computes the deviance -- with gamma functions -- for q-1 parameters
##############################################################################

dm_devG <- function(prop, disp, y){
  ## prop has length of q-1
  ## disp has length 1
  ## y has q rows and n columns 
  
  prop <- c(prop, 1 - sum(prop))
  
  ll_mod <- sum(lgamma(y + prop * disp) - lgamma(prop * disp) , na.rm = TRUE )
  
  prop_sat <- y/matrix(colSums(y), nrow(y), ncol(y), byrow = TRUE)
  
  ll_sat <- sum(lgamma(y + prop_sat * disp) - lgamma(prop_sat * disp) , 
    na.rm = TRUE) # Inf for y = 0
  
  D <- 2 * (ll_sat - ll_mod)
  
  return(D)
  
}
