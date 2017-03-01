
dm_CRadjustmentOneGroup <- function(y, disp, prop){
  # y martix q x n
  # prop vector of length q
  # If something is wrong, return NAs
  
  q <- nrow(y)
  
  # Check for 0s in rows (features)
  keep_row <- rowSums(y) > 0
  
  # Last feature can not be zero since
  # we use the last feature as a denominator in logit
  if(keep_row[q] == 0)
    return(NA)
  
  y <- y[keep_row, , drop=FALSE]
  q <- nrow(y)
  
  # Check for 0s in cols (replicates)
  keep_col <- colSums(y) > 0
  y <- y[, keep_col, drop=FALSE]
  
  n <- ncol(y) 
  prop <- prop[keep_row]
  
  adj <- log(det(n * (- dm_HessianG(prop = prop[-q], disp, y)) ))/2 
  ## with Gamma functions ## if prop is NULL then:
  # Error in is.data.frame(x) :
  # dims [product 6] do not match the length of object [0]
  
  return(adj)
  
}
