#' @importFrom ggplot2 ggplot aes_string theme_bw xlab ylab geom_histogram theme
#'   element_text coord_cartesian geom_text geom_bar scale_fill_manual


dm_plotDataFeatures <- function(tt){
  
  df <- data.frame(tt = tt)

  ggp <- ggplot(df, aes_string(x = "tt")) +
    theme_bw() +
    xlab("Number of features per gene") +
    ylab("Frequency") +
    geom_histogram(fill = "seagreen4", 
      breaks = seq(min(df$tt), max(df$tt), by = 1)) +
    theme(axis.text = element_text(size=16), 
      axis.title = element_text(size=18, face="bold"), 
      plot.title = element_text(size=18, face="bold")) +
    coord_cartesian(xlim = c(0, max(tt) + 2)) +
    geom_text(data = data.frame(x = Inf, y = Inf, 
      label = paste0(length(tt), " genes   \n ", sum(tt) , " features   ")), 
      aes_string(x = "x", y = "y", label = "label"), hjust = 1, vjust = 2, 
      size = 6)
  
  
  return(ggp)
  
}




dm_plotDataBlocks <- function(tt){
  
  
  df <- data.frame(tt = tt)
  binwidth <- ceiling(max(df$tt)/100)
  
  ggp <- ggplot(df, aes_string(x = "tt")) +
    theme_bw() +
    xlab("Number of blocks per gene") +
    ylab("Frequency") +
    geom_histogram(fill = "mediumpurple4", 
      breaks = seq(min(df$tt), max(df$tt), by = binwidth)) +
    theme(axis.text = element_text(size=16), 
      axis.title = element_text(size=18, face="bold"), 
      plot.title = element_text(size=18, face="bold")) +
    coord_cartesian(xlim = c(0, max(tt) + 2)) +
    geom_text(data = data.frame(x = Inf, y = Inf, 
      label = paste0(length(tt), " genes   \n", sum(tt) , " blocks   ")), 
      aes_string(x = "x", y = "y", label = "label"), hjust = 1, vjust = 2, 
      size = 6)
  
  return(ggp)
  
  
}



dm_plotDataSnps <- function(tt){
  
  df <- data.frame(tt = tt)
  binwidth <- ceiling(max(df$tt)/100)
  
  ggp <- ggplot(df, aes_string(x = "tt")) +
    theme_bw() +
    xlab("Number of SNPs per gene") +
    ylab("Frequency") +
    geom_histogram(fill = "royalblue4", 
      breaks = seq(min(df$tt), max(df$tt), by = binwidth)) +
    theme(axis.text = element_text(size=16), 
      axis.title = element_text(size=18, face="bold"), 
      plot.title = element_text(size=18, face="bold")) +
    coord_cartesian(xlim = c(0, max(tt) + 2)) +
    geom_text(data = data.frame(x = Inf, y = Inf, 
      label = paste0(length(tt), " genes   \n", sum(tt) , " SNPs   ")), 
      aes_string(x = "x", y = "y", label = "label"), hjust = 1, vjust = 2, 
      size = 6)
  
  
  return(ggp)
  
  
}






#' Plot the frequency of present features
#' 
#' @param info Data frame with \code{gene_id} and \code{feature_id} of ALL
#'   features
#' @param ds_info Data frame with \code{gene_id} and \code{feature_id} of ONLY
#'   DS features
#'   
#' @return \code{ggplot} object

dm_plotDataDSInfo <- function(info, ds_info){
  
  ds_info_spl <- split(as.character(ds_info$feature_id), factor(ds_info$gene_id))
  info_spl <- split(as.character(info$feature_id), factor(info$gene_id))
  
  genes <- names(info_spl)
  
  tas <- table(unlist(lapply(names(ds_info_spl), function(g){
    # g = "FBgn0004636"
    
    if(! g %in% genes)
      return(NA)
    
    features_in <- sum(info_spl[[g]] %in% ds_info_spl[[g]])
    return(features_in)
    
  })), useNA = "always")
  
  tas <- tas[c(length(tas), 1:(length(tas) -1))]
  names(tas)[1] <- "NoGene"
  
  
  df <- data.frame(x = factor(names(tas), levels = names(tas)), 
    y = as.numeric(tas), colors = "2", stringsAsFactors = FALSE)
  df[df$x == "NoGene" | df$x == "0", "colors"] <- "1"
  df$colors <- factor(df$colors)
  
  ggp <- ggplot(data = df, 
    aes_string(x = "x", y = "y", label = "y", fill = "colors")) +
    geom_bar(stat = "identity") +
    geom_text(hjust = 0.5, vjust = 0, size = 6) +
    theme_bw() +
    xlab("Number of DS features left within DS gene") +
    ylab("Number of DS genes") +
    theme(axis.text = element_text(size=16), 
      axis.title = element_text(size=18, face="bold"), 
      plot.title = element_text(size=18, face="bold"), 
      legend.position = "none") +
    scale_fill_manual(values = c("darkred", "grey"))
  
  
  return(ggp)
  
}





