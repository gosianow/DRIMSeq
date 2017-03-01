##############################################################################
# Score-function for q-1 parameters -- with gamma functions
##############################################################################

## Using gamma functions
dm_scoreG <- function(prop, disp, y){  
  
  q <- nrow(y)
  N <- ncol(y)
  yq <- y[q,]
  y <- y[-q, , drop=FALSE]
  propq <- 1-sum(prop)
  
  S <- disp * rowSums( digamma(y + prop * disp) - 
      digamma(prop * disp) - matrix(digamma(yq + disp * propq) - 
          digamma(disp * propq), nrow = q-1, ncol = N, byrow = TRUE) ) 
  
  return(S)
  
} 

dm_scoreG_neg <- function(prop, disp, y)
  -dm_scoreG(prop, disp, y)  


dm_score_regG <- function(b, x, disp, y){  
  
  y <- t(y)
  
  n <- nrow(x)
  q <- ncol(y)
  p <- ncol(x)
  
  b <- matrix(b, p, q-1)
  
  z <- exp(x %*% b)
  prop <- z/(1 + rowSums(z))
  propq <- 1 - rowSums(prop)
  
  yq <- y[, q]
  y <- y[, -q, drop = FALSE]
  
  
  S <- t(x) %*% ((digamma(y + prop*disp) - 
      digamma(prop*disp) + digamma(yq + propq*disp) - 
      digamma(propq*disp)) * prop * (1 - prop) * disp)
  
  
  return(as.numeric(S))
  
  
  
} 

dm_score_regG_neg <- function(b, design, disp, y)
  -dm_score_regG(b, x = design, disp, y)

##############################################################################
# Score-function for q-1 parameters -- with sums
##############################################################################

dm_score <- function(prop, disp, y){  
  
  q <- nrow(y)
  N <- ncol(y) 
  yq <- y[q, ]
  y <- y[-q, , drop=FALSE]
  propq <- 1-sum(prop)
  
  S <- rep(0, q-1)
  
  for(j in 1:N){ 
    # j=1
    if(yq[j] == 0) Sqj <- 0
    else Sqj <- sum(disp / (propq * disp + 1:yq[j] - 1))    
    
    for(i in 1:(q-1)){
      # i=1
      if(y[i,j] == 0) Sij <- 0
      else Sij <- sum(disp / (prop[i] * disp + 1:y[i,j] - 1)) 
      S[i] <- S[i] + Sij
    }
    S <- S - Sqj
  }
  
  return(S)
  
}











