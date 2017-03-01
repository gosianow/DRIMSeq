##############################################################################
## Computes the DM log-likelihood -- with gamma functions -- for q-1 parameters
##############################################################################

dm_likG <- function(prop, prec, y){
  # prop has length of q-1
  # prec has length 1
  # y has q rows and n columns
  # This function returns likelihhod without normalizing component, 
  # but it is OK for optimization and the LR test
  
  n <- ncol(y)
  m <- colSums(y)
  
  prop <- c(prop, 1 - sum(prop))
  
  l <- n * lgamma(prec) - sum(lgamma(m + prec)) + 
    sum( colSums( lgamma(y + prop * prec) - lgamma(prop * prec) ) )
  
  # normalizing_part <- sum(lgamma(m + 1) - colSums(lgamma(y + 1)))
  
  return(l)
  
}

dm_likG_neg <- function(prop, prec, y) 
  -dm_likG(prop, prec, y)



dm_lik_regG_prop <- function(y, prec, prop){
  # y n x q matrix !!!
  # prop n x q matrix of fitted proportions
  
  n <- nrow(y)
  
  m <- rowSums(y)
  
  l <- n * lgamma(prec) - sum(lgamma(m + prec)) + 
    sum(rowSums(lgamma(y + prec * prop) - lgamma(prop * prec)))
  
  # normalizing_part <- sum(lgamma(m + 1) - rowSums(lgamma(y + 1)))
  
  return(l)
  
}


dm_lik_regG <- function(b, x, prec, y){
  ## b has length of (q-1) * p
  ## x is a matrix n x p
  ## prec has length 1
  ## y q x n matrix
  ## This function returns likelihhod without normalizing component, 
  ## but it is OK for optimization and the LR test
  
  y <- t(y) # n x q
  
  # Get prop from x and b
  q <- ncol(y)
  p <- ncol(x)
  
  b <- matrix(b, p, q-1)
  
  z <- exp(x %*% b)
  prop <- z/(1 + rowSums(z)) # n x (q-1)
  prop <- cbind(prop, 1 - rowSums(prop))
  
  l <- dm_lik_regG_prop(y = y, prec = prec, prop = prop)

  return(l)
  
}

dm_lik_regG_neg <- function(b, design, prec, y) 
  -dm_lik_regG(b, x = design, prec, y)

##############################################################################
## Computes the DM log-likelihood -- with sums -- for q-1 parameters
##############################################################################

dm_lik <- function(prop, prec, y){
  # prop has length q-1
  # prec has length 1
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
    
    l <- l - sum(log(prec + 1:m[i] - 1))  
    
    for(j in 1:q){   
      # j=3
      
      if(y[j, i] == 0){
        lji <- 0
      }else{
        lji <- sum(log(prop[j] * prec + 1:y[j, i] - 1)) 
      }     
      
      l <- l + lji      
      
    }
    
  }
  
  return(l)
  
}


##############################################################################
## Computes the BB log-likelihood -- with gamma functions -- for q parameters
##############################################################################


bb_likG <- function(prop, prec, y){
  # prop has length of q
  # prec has length 1
  # y has q rows and n number of columns
  
  m <- colSums(y)
  q <- length(prop)
  
  l <- rep(NA, q)
  
  for(i in 1:q){
    
    l[i] <- dm_likG(prop = prop[i], prec = prec, 
      y = rbind(y[i, ], m - y[i, ]))
    
  }
  
  # l vector of length q
  return(l)
  
}


  
bb_lik_regG_prop <- function(y, prec, prop){
  # y n x q matrix !!!
  # prop n x q matrix of fitted proportions
  
  m <- rowSums(y)
  q <- ncol(prop)
  
  l <- rep(NA, q)
  
  for(i in 1:q){
    
    l[i] <- dm_lik_regG_prop(y = cbind(y[, i], m - y[, i]), prec = prec, 
      prop = cbind(prop[, i], 1 - prop[, i]))
    
  }
  
  # l vector of length q
  return(l)
  
}





