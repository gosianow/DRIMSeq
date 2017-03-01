#' @importFrom stats rgamma

dm_rdirichlet <- function(n = 1, alpha){
  Gam <- matrix(0,n,length(alpha))
  for(i in 1:length(alpha)) Gam[,i] <- rgamma(n, shape=alpha[i])
  Gam/rowSums(Gam)
}