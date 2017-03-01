##############################################################################
# Hessian for q-1 parameters -- with gamma functions
##############################################################################

dm_HessianG <- function(prop, disp, y){  
  # prop has length of q-1
  # disp has length 1
  # y has q rows and n columns
  
  q <- nrow(y)
  n <- ncol(y)
  Djj <- rep(0, q-1)
  propq <- 1-sum(prop)
  yq <- y[q, ]
  y <- y[-q, , drop=FALSE]
  
  Djl <- disp^2 * sum(trigamma(yq + disp * propq) - trigamma(disp * propq))
  
  Djj <- disp^2 * rowSums(trigamma(y + prop * disp) - trigamma(prop * disp))
  
  H <- matrix(Djl, q-1, q-1)
  
  diag(H) <- diag(H) + Djj
  
  # H martix (q-1) x (q-1)
  return(H)
  
}






##############################################################################
# Hessian for q-1 parameters -- with sums
##############################################################################

dm_Hessian <- function(prop, disp, y){  
  # prop has length of q-1
  # disp has length 1
  # y has q rows and n columns
  
  q <- nrow(y)
  n <- ncol(y)
  propq <- 1 - sum(prop)
  
  Djj <- rep(0, q-1)
  Djl <- 0
  
  for(i in 1:n){
    # i=1
    
    if(y[q, i] == 0){
      Djl <- Djl + 0
    }else{
      Djl <- Djl +  sum(-disp^2 / (propq * disp + 1:y[q, i] - 1) ^2) 
    } 
    
    for(j in 1:(q-1)){
      # j=1
      
      if(y[j,i] == 0){
        Djj[j] <- Djj[j] + 0
      }else{
        Djj[j] <- Djj[j] + sum(-disp^2 / (prop[j] * disp + 1:y[j,i] - 1) ^2)
      }  
      
    }
  }
  
  H <- matrix(Djl, q-1, q-1)
  
  diag(H) <- diag(H) + Djj
  
  # H martix (q-1) x (q-1)
  return(H)
  
}




