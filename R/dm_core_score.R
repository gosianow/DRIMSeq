##############################################################################
# Score-function for q-1 parameters -- with gamma functions
##############################################################################


dm_scoreG <- function(prop, disp, y){  
  # prop has length of q-1
  # disp has length 1
  # y has q rows and n columns
  
  q <- nrow(y)
  n <- ncol(y)
  yq <- y[q,]
  y <- y[-q, , drop=FALSE]
  propq <- 1-sum(prop)
  
  S <- disp * rowSums( digamma(y + prop * disp) - 
      digamma(prop * disp) - matrix(digamma(yq + disp * propq) - 
          digamma(disp * propq), nrow = q-1, ncol = n, byrow = TRUE) ) 
  
  return(S)
  
} 

dm_scoreG_neg <- function(prop, disp, y)
  -dm_scoreG(prop, disp, y)  


dm_score_regG <- function(b, x, disp, y){  
  ## b has length of (q-1) * p
  ## x is a matrix n x p
  ## disp has length 1
  ## y has q rows and n columns
  
  y <- t(y) # n x q
  
  n <- nrow(x)
  q <- ncol(y)
  p <- ncol(x)
  
  b <- matrix(b, p, q-1) # p x (q-1)
  
  z <- exp(x %*% b)
  prop_qm1 <- z/(1 + rowSums(z))
  prop <- cbind(prop_qm1, 1 - rowSums(prop_qm1))
  
  y_qm1 <- y[, -q, drop = FALSE] # n x (q-1)
  
  S <- t(x) %*% 
    (disp * prop_qm1 * (- rowSums(digamma(y + prop*disp) * prop) +
        digamma(y_qm1 + prop_qm1*disp) +
        rowSums(digamma(prop*disp) * prop) -
        digamma(prop_qm1*disp))) # p x (q-1)
  
  return(c(S))
  
} 


dm_score_regG_neg <- function(b, design, disp, y)
  -dm_score_regG(b, x = design, disp, y)


##############################################################################
# Score-function for q-1 parameters -- with sums
##############################################################################

dm_score <- function(prop, disp, y){  
  # prop has length of q-1
  # disp has length 1
  # y has q rows and n columns
  
  q <- nrow(y)
  n <- ncol(y) 
  propq <- 1 - sum(prop)
  
  S <- rep(0, q-1)
  
  for(i in 1:n){ 
    # i=1
    
    if(y[q, i] == 0){
      Sqi <- 0
    }else{
      Sqi <- sum(disp / (propq * disp + 1:y[q, i] - 1))
    }     
    
    for(j in 1:(q-1)){
      # j=1
      
      if(y[j, i] == 0){
        Sji <- 0
      }else{
        Sji <- sum(disp / (prop[j] * disp + 1:y[j, i] - 1))
      }  
      
      S[j] <- S[j] + Sji
      
    }
    S <- S - Sqi
  }
  
  return(S)
  
}











