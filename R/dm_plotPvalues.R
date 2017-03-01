#' @importFrom ggplot2 ggplot aes_string theme_bw xlab ylab geom_histogram theme
#'   element_text coord_cartesian geom_text

dm_plotPValues <- function(pvalues){
  
  df <- data.frame(pvalues = pvalues[!is.na(pvalues)])
  
  ggp <- ggplot(df, aes_string(x = "pvalues")) +
    theme_bw() +
    xlab("P-Values") +
    ylab("Frequency") +
    geom_histogram(breaks = seq(0, 1, by = 0.01), fill = "deeppink4") +
    theme(axis.text = element_text(size=16), 
      axis.title = element_text(size=18, face="bold"), 
      plot.title = element_text(size=16, face="bold")) +
    coord_cartesian(xlim = c(0, 1)) +
    geom_text(data = data.frame(x = Inf, y = Inf, 
      label = paste0(nrow(df), " tests       ")), 
      aes_string(x = "x", y = "y", label = "label"), 
      hjust = 1, vjust = 3, size = 6)
  
  return(ggp)
  
}










