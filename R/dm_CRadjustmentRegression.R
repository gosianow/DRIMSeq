

dm_CRadjustmentRegression <- function(y, x, prec, prop){  
  # prop q x n matrix of fitted proportions
  # y q x n matrix
  # x n x p matrix with the design
  # y can not have any rowSums(y) == 0 - assured during dmFilter

  if(any(is.na(prop[1, ])))
    return(NA)
  
  prop <- t(prop) # n x q
  y <- t(y) # n x q
  
  n <- nrow(y)
  
	H <- dm_Hessian_regG_prop(y = y, prec = prec, prop = prop, x = x)

  adj <- log(det(n * (-H))) / 2 
  
  if(is.na(adj))
    return(NA) 

  if(abs(adj) == Inf)
    return(NA) 
  
  return(adj)
  
}





















