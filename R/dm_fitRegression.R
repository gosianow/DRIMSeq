
#' @importFrom stats optim nlminb

dm_fitRegression <- function(y, design, 
  disp, coef_mode = "optim", coef_tol = 1e-12){
  # y can not have any rowSums(y) == 0 - assured during dmFilter
  
  q <- nrow(y)
  p <- ncol(design)
  n <- ncol(y)
  
  # NAs for genes with one feature
  if(q < 2 || is.na(disp)){
    b <- matrix(NA, nrow = q, ncol = p)
    prop <- matrix(NA, nrow = q, ncol = n)
    rownames(prop) <- rownames(y)
    rownames(b) <- rownames(y)
    return(list(b = b, lik = NA, fit = prop))  
  }
  
  # Get the initial values for b 
  # Add 1 to get rid of NaNs and Inf and -Inf values in logit
  # Double check if this approach is correct!!!
  # Alternative, use 0s: b_init = rep(0, p*(q-1))
  yt <- t(y) + 1 # n x q
  prop <- yt / rowSums(yt) # n x q
  logit_prop <- log(prop / prop[, q])
  b_init <- c(MASS::ginv(design) %*% logit_prop[, -q, drop = FALSE])
  
  
  switch(coef_mode, 
    
    optim = { 
      
      # Minimization
      co <- optim(par = b_init, fn = dm_lik_regG_neg, gr = dm_score_regG_neg,
        design = design, disp = disp, y = y,
        method = "BFGS",
        control = list(reltol = coef_tol))
      
      if(co$convergence == 0){
        b <- rbind(t(matrix(co$par, p, q-1)), rep(0, p))
        lik <- -co$value
      }else{
        b <- matrix(NA, nrow = q, ncol = p)
        lik <- NA
      }
      
    },
    
    nlminb = { 
      
      # Minimization
      co <- nlminb(start = b_init, objective = dm_lik_regG_neg, 
        gradient = dm_score_regG_neg, hessian = NULL,
        design = design, disp = disp, y = y,
        control = list(rel.tol = coef_tol))
      
      if(co$convergence == 0){
        b <- rbind(t(matrix(co$par, p, q-1)), rep(0, p))
        lik <- -co$objective
      }else{
        b <- matrix(NA, nrow = q, ncol = p)
        lik <- NA
      }
      
    },
    
    Rcgmin = {
      
      # Minimization
      # Can not use x as an argument because grad() uses x
      co <- Rcgmin::Rcgmin(par = b_init, fn = dm_lik_regG_neg, 
        gr = dm_score_regG_neg, 
        design = design, disp = disp, y = y) 
      
      if(co$convergence == 0){
        b <- rbind(t(matrix(co$par, p, q-1)), rep(0, p))
        lik <- -co$value
      }else{
        b <- matrix(NA, nrow = q, ncol = p)
        lik <- NA
      }
      
    })
  
  # Compute the fitted proportions
  if(!is.na(lik)){
    z <- exp(design %*% t(b[-q, , drop = FALSE]))
    prop <- z/(1 + rowSums(z)) # n x (q-1)
    prop <- t(cbind(prop, 1 - rowSums(prop))) # q x n
  }else{
    prop <- matrix(NA, nrow = q, ncol = n)
  }
  
  rownames(prop) <- rownames(y)
  rownames(b) <- rownames(y)
  
  # b matrix q x p
  # lik vector of length 1
  # fit matrix q x n
  return(list(b = b, lik = lik, fit = prop))  
  
}


# -----------------------------------------------------------------------------
# Fitting the Beta-binomial model
# Currently, recalculating the BB likelihoods and coefficients using the 
# DM fittings/proportions


bb_fitRegression <- function(y, design, disp, fit){
  # y can not have any rowSums(y) == 0 - assured during dmFilter
  
  q <- nrow(y)
  p <- ncol(design)
  n <- ncol(y)
  
  # NAs for genes with one feature
  if(q < 2 || is.na(disp)){
    b <- matrix(NA, nrow = q, ncol = p)
    prop <- matrix(NA, nrow = q, ncol = n)
    rownames(prop) <- rownames(y)
    rownames(b) <- rownames(y)
    return(list(b = b, lik = rep(NA, q), fit = prop))  
  }
  
  y <- t(y) # n x q
  prop <- t(fit) # n x q
  
  lik <- bb_lik_regG_prop(y = y, disp = disp, prop = prop)

  # Get the coefficients like in edgeR::mglmOneWay
  # But use MASS::ginv instead of solve since the design does not have to be 
  # a squared matrix
  
  # Keep only the unique rows in the design and logit_prop
  unique_samps <- !duplicated(design)
  
  design <- design[unique_samps, , drop = FALSE]
  prop <- prop[unique_samps, , drop = FALSE]
  
  logit_prop <- log(prop / (1 - prop))

  b <- t(MASS::ginv(design) %*% logit_prop)

  rownames(b) <- rownames(y)
  
  # b matrix q x p
  # lik vector of length q
  # fit matrix q x n
  return(list(b = b, lik = lik, fit = fit))  
  
}






























