##############################################################################
# Hessian for Dirichlet-multinomial model gamma+ and k-1 proportions 
##############################################################################

# prop <- prop[-length(prop)]

### Using gamma functions
dm_HessianG <- function(prop, disp, y){  
  ## prop has length of k-1
  
  k <- nrow(y)
  N <- ncol(y)
  D <- rep(0, k-1)
  Dil <- 0
  propk <- 1-sum(prop)
  ykm1 <- y[-k, , drop=FALSE]
  
  Dil <- disp^2 * sum(trigamma(y[k,] + disp * propk) - trigamma(disp * propk))
  
  D <- disp^2 * rowSums(trigamma(ykm1 + prop * disp) - trigamma(prop * disp))
  
  H <- matrix(Dil, k-1, k-1)
  diag(H) <- diag(H) + D
  
  return(H)
  
  
}


dm_Hessian <- function(prop, disp, y){  
  ## prop has length of k-1
  
  k <- nrow(y)
  N <- ncol(y)
  D <- rep(0, k-1)
  Dil <- 0
  propk <- 1-sum(prop)
  
  for(j in 1:N){
    # j=1
    if(y[k, j] == 0) Dil <- Dil + 0
    else Dil <- Dil +  sum(-disp^2 / (propk * disp + 1:y[k, j] - 1) ^2) 
    
    for(i in 1:(k-1)){
      # i=1
      if(y[i,j] == 0) Dii <- 0
      else Dii <- sum(-disp^2 / (prop[i] * disp + 1:y[i,j] - 1) ^2) 
      D[i] <- D[i] + Dii
    }
  }
  
  H <- matrix(Dil, k-1, k-1)
  diag(H) <- diag(H) + D
  
  return(H)
  
}




