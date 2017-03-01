##############################################################################
## Computes the log-likelihood -- with gamma functions -- for q-1 parameters
##############################################################################

dm_likG <- function(prop, disp, y){
  # prop has length of q-1
  # disp has length 1
  # y has q rows and n columns
  # This function returns likelihhod without normalizing component, 
  # but it is OK for optimization and the LR test
  
  n <- ncol(y)
  m <- colSums(y)
  
  prop <- c(prop, 1 - sum(prop))
  
  l <- n * lgamma(disp) - sum(lgamma(m + disp)) + 
    sum( colSums( lgamma(y + prop * disp) - lgamma(prop * disp) ) )
  
  # normalizing_part <- sum(lgamma(m + 1) - colSums(lgamma(y + 1)))
  
  return(l)
  
}

dm_likG_neg <- function(prop, disp, y) 
  -dm_likG(prop, disp, y)


bb_likG <- function(prop, disp, y){
  # prop has length of q
  # disp has length 1
  # y has q rows and n number of columns
  
  m <- colSums(y)
  
  l <- rep(NA, length(prop))
  
  for(i in 1:length(prop)){
    
    l[i] <- dm_likG(prop[i], disp, y = rbind(y[i, ], m - y[i, ]))
    
  }
  
  return(l)
  
}




dm_lik_regG <- function(b, x, disp, y){
  ## b has length of (q-1) * p
  ## x is a matrix n x p
  ## disp has length 1
  ## y has q rows and n columns
  ## This function returns likelihhod without normalizing component, 
  ## but it is OK for optimization and the LR test
  
  y <- t(y) # n x q
  n <- nrow(x)
  
  # Get prop from x and b
  q <- ncol(y)
  p <- ncol(x)
  
  b <- matrix(b, p, q-1)
  
  z <- exp(x %*% b)
  prop <- z/(1 + rowSums(z)) # n x (q-1)
  prop <- cbind(prop, 1 - rowSums(prop))
  
  m <- rowSums(y)
  
  l <- sum(rowSums(lgamma(y + disp * prop) - lgamma(prop * disp)))
  
  l <- n * lgamma(disp) - sum(lgamma(m + disp)) + l
  
  # normalizing_part <- sum(lgamma(m + 1) - rowSums(lgamma(y + 1)))
  
  return(l)
  
}

dm_lik_regG_neg <- function(b, design, disp, y) 
  -dm_lik_regG(b, x = design, disp, y)

##############################################################################
## Computes the log-likelihood -- with sums -- for q-1 parameters
##############################################################################

dm_lik <- function(prop, disp, y){
  # prop has length q-1
  # disp has length 1
  # y has q rows and n columns
  # This function returns likelihhod without normalizing component, 
  # but it is OK for optimization and the LR test
  
  q <- nrow(y)
  n <- ncol(y)
  m <- colSums(y) 
  prop <- c(prop, 1 - sum(prop)) 
  
  l <- 0
  
  for(i in 1:n){  
    # i=1
    
    l <- l - sum(log(disp + 1:m[i] - 1))  
    
    for(j in 1:q){   
      # j=3
      
      if(y[j, i] == 0){
        lji <- 0
      }else{
        lji <- sum(log(prop[j] * disp + 1:y[j, i] - 1)) 
      }     
      
      l <- l + lji      
      
    }
    
  }
  
  return(l)
  
}

