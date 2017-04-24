##############################################################################
# Hessian for q-1 parameters -- with gamma functions
##############################################################################

dm_HessianG <- function(prop, prec, y){  
  # prop has length of q-1
  # prec has length 1
  # y has q rows and n columns
  
  q <- nrow(y)
  n <- ncol(y)
  Djj <- rep(0, q-1)
  propq <- 1-sum(prop)
  yq <- y[q, ]
  y <- y[-q, , drop=FALSE]
  
  Djl <- prec^2 * sum(trigamma(yq + propq * prec) - trigamma(propq * prec))
  
  Djj <- prec^2 * rowSums(trigamma(y + prop * prec) - trigamma(prop * prec))
  
  H <- matrix(Djl, q-1, q-1)
  
  diag(H) <- diag(H) + Djj
  
  # H martix (q-1) x (q-1)
  return(H)
  
}



dm_Hessian_regG_prop <- function(y, prec, prop, x){
  # y n x q matrix !!!
  # prop n x q matrix of fitted proportions
  # x n x p matrix with the design
  
  q <- ncol(y)
  p <- ncol(x)
  
  prop_prec <- prop * prec
  
  # Some calculations for later
  ldg_y <- digamma(y + prop_prec)
  ldg <- digamma(prop_prec)
  
  ldgp_y <- ldg_y * prop
  ldgp_y_sum <- rowSums(ldgp_y)
  ldgp <- ldg * prop
  ldgp_sum <- rowSums(ldgp)
  
  ltgp_y <- trigamma(y + prop_prec) * prop * prec
  ltgp <- trigamma(prop_prec) * prop * prec
  
  ltgpp_y <- ltgp_y * prop
  ltgpp_y_sum <- rowSums(ltgpp_y)
  ltgpp <- ltgp * prop
  ltgpp_sum <- rowSums(ltgpp)
  
  H <- matrix(0, p*(q-1), p*(q-1))
  
  rownames(H) <- colnames(H) <- paste0(rep(colnames(y)[-q], each = p), 
    ":", rep(colnames(x), q-1))
  
  jp_index <- matrix(1:p, nrow = q-1, ncol = p, byrow = TRUE) + p * 0:(q-2)
  
  for(jp in 1:(q-1)){
    
    for(jpp in 1:jp){
      # jp = 1; jpp = 1
      
      W <- prec * ( - prop[, jp] * prop[, jpp] * 
          (-ldgp_y_sum + ldg_y[, jp] + ldgp_sum - ldg[, jp]) +
          prop[, jp] * 
          (ltgpp_y_sum * prop[, jpp] - ltgpp_y[, jpp] +
              ldgp_y_sum * prop[, jpp] - ldgp_y[, jpp] -
              ltgp_y[, jp] * prop[, jpp] -
              ltgpp_sum * prop[, jpp] + ltgpp[, jpp] -
              ldgp_sum * prop[, jpp] + ldgp[, jpp] +
              ltgp[, jp] * prop[, jpp]) )
      
      if(jp == jpp){
        
        W <- W + prec * ( prop[, jp] *
            (-ldgp_y_sum + ldg_y[, jp] + ldgp_sum - ldg[, jp]) +
            prop[, jp] * (ltgp_y[, jp] + ltgp[, jp]) )
        
      }
      
      h <- t(x) %*% (W * x)
      
      H[jp_index[jp, ], jp_index[jpp, ]] <- h
      
      # Use the fact that H is symetric
      if(!jp == jpp){
        H[jp_index[jpp, ], jp_index[jp, ]] <- t(h)
      }
      
    }
    
  }
  
  # H martix p(q-1) x p(q-1)
  return(H)
  
}


dm_Hessian_regG <- function(b, x, prec, y){  
  ## b has length of (q-1) * p
  ## x is a matrix n x p
  ## y has q rows and n columns
  
  y <- t(y) # n x q
  
  q <- ncol(y)
  p <- ncol(x)
  
  b <- matrix(b, p, q-1) # p x (q-1)
  
  z <- exp(x %*% b) # n x (q-1)
  
  prop_qm1 <- z/(1 + rowSums(z))
  
  prop <- cbind(prop_qm1, 1 - rowSums(prop_qm1)) # n x q
  
  H <- dm_Hessian_regG_prop(y = y, prec = prec, prop = prop, x = x)
  
  return(H)
  
} 



# Hessian for the multinomial distribution 

m_Hessian_regG <- function(b, x, y){  
  ## b has length of (q-1) * p
  ## x is a matrix n x p
  ## y has q rows and n columns
  
  y <- t(y) # n x q
  
  q <- ncol(y)
  p <- ncol(x)
  
  b <- matrix(b, p, q-1) # p x (q-1)
  
  z <- exp(x %*% b) # n x (q-1)
  
  prop_qm1 <- z/(1 + rowSums(z))
  
  prop <- cbind(prop_qm1, 1 - rowSums(prop_qm1)) # n x q
  
  m <- rowSums(y) # n 
  
  H <- matrix(0, p*(q-1), p*(q-1))
  
  rownames(H) <- colnames(H) <- paste0(rep(colnames(y)[-q], each = p), 
    ":", rep(colnames(x), q-1))
  
  jp_index <- matrix(1:p, nrow = q-1, ncol = p, byrow = TRUE) + p * 0:(q-2)
  
  for(jp in 1:(q-1)){
    
    for(jpp in 1:jp){
      # jp = 1; jpp = 1
      
      W <- m * prop[, jp] * (prop[, jpp] - as.numeric(jp == jpp))
      
      h <- t(x) %*% (W * x)
      
      H[jp_index[jp, ], jp_index[jpp, ]] <- h
      
      # Use the fact that H is symetric
      if(!jp == jpp){
        H[jp_index[jpp, ], jp_index[jp, ]] <- t(h)
      }
      
    }
    
  }
  
  # H martix p(q-1) x p(q-1)
  return(H)
  
} 




##############################################################################
# Hessian for q-1 parameters -- with sums
##############################################################################

dm_Hessian <- function(prop, prec, y){  
  # prop has length of q-1
  # prec has length 1
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
      Djl <- Djl +  sum(-prec^2 / (propq * prec + 1:y[q, i] - 1) ^2) 
    } 
    
    for(j in 1:(q-1)){
      # j=1
      
      if(y[j,i] == 0){
        Djj[j] <- Djj[j] + 0
      }else{
        Djj[j] <- Djj[j] + sum(-prec^2 / (prop[j] * prec + 1:y[j,i] - 1) ^2)
      }  
      
    }
  }
  
  H <- matrix(Djl, q-1, q-1)
  
  diag(H) <- diag(H) + Djj
  
  # H martix (q-1) x (q-1)
  return(H)
  
}




