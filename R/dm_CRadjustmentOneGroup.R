
dm_CRadjustmentOneGroup <- function(y, prec, prop){
  # y martix q x n
  # prop vector of length q
  # If something is wrong, return NAs
  
  if(any(is.na(prop)))
    return(NA)
  
  q <- nrow(y)
  
  # NAs for genes with one feature
  if(q < 2 || is.na(prec)) 
    return(NA)
  
  # Check for 0s in rows (features)
  keep_row <- rowSums(y) > 0
  # There must be at least two non-zero features
  if(sum(keep_row) < 2) 
    return(NA)
  
  # Last feature can not be zero since
  # we use the last feature as a denominator in logit
  if(keep_row[q] == 0)
    return(NA)
  
  y <- y[keep_row, , drop=FALSE]
  q <- nrow(y)
  
  # Check for 0s in columns (replicates)
  keep_col <- colSums(y) > 0
  y <- y[, keep_col, drop=FALSE]
  
  prop <- prop[keep_row]
  n <- ncol(y) 
  
  H <- dm_HessianG(prop = prop[-q], prec, y)
  
  adj <- log(det(n * (- H) ))/2 
  ## with Gamma functions ## if prop is NULL then:
  # Error in is.data.frame(x) :
  # dims [product 6] do not match the length of object [0]
  
  return(adj)
  
}
