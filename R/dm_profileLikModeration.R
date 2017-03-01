##############################################################################
# calculate the moderated profile likelihood
##############################################################################

#' @importFrom stats loess predict loess.control

dm_profileLikModeration <- function(loglik, mean_expression, 
  disp_moderation = "trended", disp_prior_df, disp_span){
  
  disp_grid_length <- ncol(loglik)
  
  ### Check where the grid is maximized 
  grid_max <- apply(loglik, 1, which.max)
  
  # In the calculation of moderation, do not take into account genes 
  # that have dispersion on the top and bottom boundry of the grid 
  # (skipp 4 last grid points and 1 first grid point)
  not_boundry <- grid_max < (disp_grid_length - 3) & grid_max > 1
  boundry_last <- grid_max == disp_grid_length
  
  ### Calculate the span of the boundry loglikelihoods
  if(sum(boundry_last) > 1){
    loglik_span_boundry <- apply(loglik[boundry_last, , drop = FALSE], 1, 
      function(x){max(x) - min(x)})
  }
  
  
  switch(disp_moderation, 
    
    common={
      
      ### Calculate the moderating likelihood
      if(sum(not_boundry) == length(not_boundry)){
        moderation <- colMeans(loglik)
      }else{
        moderation <- colMeans(loglik[not_boundry, , drop = FALSE])
      }
      
      
      # Estimate priorN - calculate the ratio between moderation lik span 
      # and lik span of boundry genes
      if(sum(boundry_last) > 10){
        
        moderation_span <- max(moderation) - min(moderation)
        span_ratio <- moderation_span / loglik_span_boundry
        priorN <- 1/span_ratio
        
        ### Use median 
        priorN <- quantile(priorN, 0.5)
        
      }else{
        priorN <- disp_prior_df
      }
      

        message(paste0("! Using ", round(priorN, 4), 
          " as a shrinkage factor !\n"))
        
        loglik <- sweep(loglik, 2, priorN * moderation, FUN = "+")

      
    },
    
    trended={
      
      moderation <- dm_movingAverageByCol(loglik = loglik, 
        mean_expression = mean_expression, not_boundry = not_boundry, 
        disp_span = disp_span)
      
      
      # Estimate priorN - calculate the ratio between moderation lik span 
      # and lik span of boundry genes
      if(sum(boundry_last) > 10){
        
        moderation_span_boundry <- apply(
          moderation[boundry_last, , drop = FALSE], 1, 
          function(x){max(x) - min(x)})
        span_ratio <- moderation_span_boundry / loglik_span_boundry
        priorN <- 1/span_ratio
        
        
        ### Do loess fitting if there is enough points. Otherwise, use median
        if(length(loglik_span_boundry) > 100){

          df_priorN_loglog <- data.frame(priorN = log10(priorN), 
            mean_expression = log10(mean_expression[boundry_last]))
          
          priorN_loess_loglog <- loess(priorN ~ mean_expression, 
            df_priorN_loglog, control = loess.control(surface = "direct"))
          priorN_predict_loglog <- predict(priorN_loess_loglog, 
            data.frame(mean_expression = log10(mean_expression)), se = FALSE)
          
          priorN <- 10 ^ priorN_predict_loglog
          
        }else{

          priorN <- quantile(priorN, 0.5)
          
        }
        
      }else{
        priorN <- disp_prior_df
      }
      
      if(length(priorN) == 1){
        message(paste0("! Using ", round(priorN, 6), 
          " as a shrinkage factor !\n"))
      }else{
        message(paste0("! Using loess fit as a shrinkage factor !\n"))
      }
      
      loglik <- loglik + priorN * moderation
      
    }
  )
  
  return(loglik)
  
}



dm_movingAverageByCol <- function(loglik, mean_expression, 
  not_boundry, disp_span){
  
  if(sum(not_boundry) == length(not_boundry)){
    
    o <- order(mean_expression)
    oo <- order(o)
    width <- floor(disp_span * nrow(loglik))
    
    moderation <- edgeR::movingAverageByCol(loglik[o,], width = width)[oo,]
    
  }else{
    
    ### Use non boundry genes for calculating the moderation
    mean_expression_not_boundry <- mean_expression[not_boundry]
    loglik_not_boundry <- loglik[not_boundry, , drop = FALSE]
    
    o <- order(mean_expression_not_boundry)
    oo <- order(o)
    
    width <- floor(disp_span * nrow(loglik_not_boundry))
    
    moderation_not_boundry <- edgeR::movingAverageByCol(
      loglik_not_boundry[o, , drop = FALSE], width = width)[oo, , drop = FALSE]
    
    ### Fill in moderation values for the boundy genes
    moderation <- matrix(NA, nrow = nrow(loglik), ncol = ncol(loglik))
    
    moderation[not_boundry, ] <- moderation_not_boundry
    
    o <- order(mean_expression)
    oo <- order(o)
    
    moderation <- moderation[o, , drop = FALSE]
    not_boundry <- not_boundry[o]
    
    ### Last value in not_boundry must be TRUE
    if(not_boundry[length(not_boundry)] == FALSE){
      
      last_true <- max(which(not_boundry))
      moderation[length(not_boundry), ] <- moderation[last_true, ]
      
      not_boundry[length(not_boundry)] <- TRUE
      
    }
    
    not_boundry_diff <- diff(not_boundry, lag = 1)
    
    not_boundry_cumsum <- cumsum(not_boundry)
    
    ### Values used for filling in the boundry NAs - swith from FALSE to TRUE
    replacement_indx <- which(not_boundry_diff == 1) + 1
    
    replaced_indx <- which(!not_boundry)
    
    replaced_freq <- as.numeric(table(not_boundry_cumsum[replaced_indx]))
    
    moderation_boundry  <- moderation[rep(replacement_indx, 
      times = replaced_freq), , drop = FALSE]
    
    moderation[!not_boundry, ] <- moderation_boundry
    
    moderation <- moderation[oo, , drop = FALSE]
    
  }
  
  
  return(moderation)
  
}









