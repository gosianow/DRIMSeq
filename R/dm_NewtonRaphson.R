dm_NewtonRaphson <- function(y, k, gamma0, verbose = FALSE){
  ### This mode is completely unprepared!!! And slow.
  
  plot = FALSE
  epsilon = 1e-05
  maxIte = 1000
  pi_initOrg <- pi_init <- rowSums(y)/sum(y)
  piMAX <- pi <- pi_init[-k]
  lik1 <- 0
  likMAX <- lik2 <- dm_lik(pi, gamma0, y)
  ite <- 1
  conv <- TRUE 
  
  # if(plot){
  #   piX <- seq(0, 1, by = 0.01)
  #   piX <- piX[c(-1, -length(piX))]
  #   loglikY <- rep(0, length(piX))
  #   for(i in 1:length(loglikY))
  #     loglikY[i] <- dm_lik(piX[i], gamma0, y)
  
  #   plot(piX, loglikY, type="l", col="deeppink", lwd=4, main = gamma0)
  #   abline(v = pi, lty =  3)
  #   plot(piX, loglikY, type="l", col="deeppink", lwd=6, main = gamma0, 
  # xlim=c(pi - 0.3, pi + 0.3))
  #   abline(v = pi, lty =  3)
  #   points(pi, lik2, pch="0")
  # }
  
  # Iterations
  while(conv & ite <= maxIte){
    if(verbose) message("\n gene:", colnames(y)[1], "gamma0:", gamma0)
    if(verbose) message("ite:", ite)
    if(verbose) message("lik2-lik1:", lik2 - lik1)
    
    if(abs(lik2 - lik1) < epsilon) conv <- FALSE
    
    score <- dm_scoreG(pi, gamma0, y)
    if(verbose) message("score:", score)
    
    # if(mode=="exp") invFIM <- dmInvExpFIMkm1(pi, gamma0, y)
    # else if(mode=="obs") invFIM <- dmInvObsFIMkm1(pi, gamma0, y)
    hessian <- dm_HessianG(pi, gamma0, y)
    
    if(verbose) message("hessian:", hessian)
    
    # Update parameter estimates
    update <- solve(hessian, score)
    pi <- pi - update
    
    
    ## check if pi is negative, then restart with new init params
    if(any(c(pi, 1-sum(pi)) <= 0 )){
      message("* Negative pi:",c(pi, 1-sum(pi)))
      ### generate new starting params
      #         randInit <- runif(k)
      #         pi_init <- randInit/sum(randInit)
      pi_init <- dm_rdirichlet(1, pi_initOrg*10)
      while(any(pi_init < 1e-10))
        pi_init <- dm_rdirichlet(1, pi_initOrg*10)
      message("pi_init:", pi_init)
      pi <- pi_init[-k]
    }
    
    if(verbose) message("pi:",c(pi, 1-sum(pi)))
    
    lik1 <- lik2
    lik2 <- dm_likG(pi, gamma0, y)
    
    if(lik2 > likMAX){
      likMAX  <- lik2
      piMAX <- pi
    }
    
    # if(plot)
    #   points(pi, lik2, pch="|")
    if(verbose) message("lik2:", lik2)
    ite <- ite + 1
  }         
  
  if(ite > maxIte){
    pi <- piMAX
    lik1 <- likMAX
  }
  
  # if(plot){
  #   points(pi, lik1, pch="*", cex=2, col="darkturquoise")
  #   dev.off()
  # }
  
  pi <- c(pi, 1-sum(pi))
  
  #   ## check if the sum of pi equals 1
  #   if(sum(pi) != 1){
  #     message("gene:", colnames(y)[1], "gamma0:", gamma0)
  #     message("**** Not 1. sum(pi) = ", sum(pi))
  #     message("**** pi:",pi)
  #   }
  
  ## check if pi is negative
  if(any(pi < 0 | pi >1 )){
    message("**** Negative pi:",pi)
  }
}

