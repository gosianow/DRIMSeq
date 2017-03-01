#' @importFrom ggplot2 ggplot aes_string theme_bw xlab ylab theme element_text
#'   guides guide_colorbar scale_colour_gradient geom_point geom_hline
#' @importFrom stats quantile na.omit


dm_plotDispersion <- function(genewise_dispersion, mean_expression, 
  nr_features = NULL, common_dispersion = NULL, low_color = "royalblue2", 
  high_color = "red2", na_value_color = "red2"){
  
  if(!is.null(nr_features)){
    
    df <- data.frame(mean_expression = log10(mean_expression + 1), 
      dispersion = log10(genewise_dispersion), nr_features = nr_features)
    
    df_quant <- min(quantile(na.omit(df$nr_features), probs = 0.95), 30)
    breaks <- seq(2, df_quant, ceiling(df_quant/10))
    
    
    ggp <- ggplot(df, aes_string(x = "mean_expression", y = "dispersion", 
      colour = "nr_features" )) +
      theme_bw() +
      xlab("Log10 of mean expression") +
      ylab("Log10 of gamma0") +
      geom_point(alpha = 0.7, na.rm = TRUE) +
      theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=18, face="bold"), 
        legend.title = element_text(size=16, face="bold"), 
        legend.text = element_text(size = 14), 
        legend.position = "top") +
      guides(colour = guide_colorbar(barwidth = 20, barheight = 0.5)) +
      scale_colour_gradient(limits = c(2, max(breaks)), 
        breaks = breaks, low = low_color, high = high_color, 
        name = "Number of features", na.value = na_value_color)
    
    
  }else{
    
    df <- data.frame(mean_expression = log10(mean_expression + 1), 
      dispersion = log10(genewise_dispersion))
    
    ggp <- ggplot(df, aes_string(x = "mean_expression", y = "dispersion")) +
      theme_bw() +
      xlab("Log10 of mean expression") +
      ylab("Log10 of gamma0") +
      geom_point(size = 1, alpha = 0.4, na.rm = TRUE) +
      theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=18, face="bold"), 
        legend.title = element_text(size=16, face="bold"), 
        legend.text = element_text(size = 14), 
        legend.position = "top")
    
  }
  
  
  if(!is.null(common_dispersion)){
    ggp <- ggp + geom_hline(yintercept = log10(common_dispersion), 
      colour = "black", linetype = "dashed")
  }
  
  
  return(ggp)
  
}

