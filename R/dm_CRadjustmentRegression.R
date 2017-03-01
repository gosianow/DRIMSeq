

dm_CRadjustmentRegression <- function(y, x, disp, prop){  
  # prop q x n matrix of fitted proportions
  # y q x n matrix
  # x n x p matrix with the design
  # y can not have any rowSums(y) == 0 - assured during dmFilter

  if(any(is.na(prop[1, ])))
    return(NA)
  
  prop <- t(prop) # n x q
  y <- t(y) # n x q
  
  adj <- log(det(n * (- dm_Hessian_regG_prop(y = y, disp = disp, 
    prop = prop, x = x)))) / 2 
  
  if(abs(adj) == Inf)
    return(NA) 
  
  return(adj)
  
}





















