
#' @importFrom stats optimHess

dm_CRadjustmentRegression <- function(y, x, prec, prop){  
  # prop q x n matrix of fitted proportions
  # y q x n matrix
  # x n x p matrix with the design
  # y can not have any rowSums(y) == 0 - assured during dmFilter
  
  if(any(is.na(prop[1, ])))
    return(NA)
  
  prop <- t(prop) # n x q
  n <- nrow(prop)
  q <- ncol(prop)
  
  H <- dm_Hessian_regG_prop(y = t(y), prec = prec, prop = prop, x = x)
  
  adj <- log(det(n * (-H))) / 2 
  
  ## If the above calculation returns NA, try optimHess()
  if(is.na(adj)){
    
    # Recalculate betas from proportions and the design
    logit_prop <- log(prop / prop[, q]) # n x q
    par <- c(MASS::ginv(x) %*% logit_prop[, -q, drop = FALSE])

    # Maximization
    H <- optimHess(par = par,
      fn = dm_lik_regG, gr = dm_score_regG,
      x = x, prec = prec, y = y)

    adj <- log(det(n * (-H))) / 2
    
    if(is.na(adj))
      return(NA) 
  }
  
  if(abs(adj) == Inf)
    return(NA) 
  
  return(adj)
  
}





















