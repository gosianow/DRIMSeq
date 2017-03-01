##############################################################################
## Computes the MoM estimate of theta (and std. error)
## similar as in dirmult package
## use as starting values for gamma0 = (1-mom)/mom
##############################################################################



dm_weirMoM <- function(y, se = FALSE){
  
  y <- y[rowSums(y) > 0, colSums(y) > 0, drop = FALSE]
  
  K <- ncol(y)
  J <- nrow(y)
  MoM <- colSums(y)/sum(y)
  Sn <- rowSums(y)
  MSP <- (J - 1)^(-1)*sum(rowSums((y/rowSums(y) - matrix(rep(MoM, J), J, K, 
    byrow = TRUE))^2)*Sn)
  MSG <- (sum(y) - J)^(-1)*sum(rowSums(y/rowSums(y)*(1 - y/rowSums(y)))*Sn)
  nc <- 1/(J - 1)*(sum(Sn) - sum(Sn^2)/sum(Sn))
  
  MoM.wh <- (MSP - MSG)/(MSP + (nc-1)*MSG)
  
  #   if(se){
  #     ## Formula by Li, ref in Weir-Hill 2002
  #     std.er <- sqrt(2*(1-MoM.wh)^2/(J-1)*((1+(nc-1)*MoM.wh)/nc)^2)
  #     list(theta = MoM.wh, se=std.er)
  #   }
  #   else MoM.wh
  
  if(MoM.wh <= 0)
    return(NA)
  
  return((1 - MoM.wh)/MoM.wh)
  
}

