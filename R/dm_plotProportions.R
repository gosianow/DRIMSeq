

dm_plotProportions_barplot <- function(prop_samp, prop_fit = NULL, 
  main = NULL, group_colors){
  
  ## Plotting
  ggp <- ggplot() +
    geom_bar(data = prop_samp, aes_string(x = "feature_id", y = "proportion", 
      group = "sample_id", fill = "group"), 
      stat = "identity", position = position_dodge(width = 0.9)) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
      axis.text=element_text(size=16), 
      axis.title = element_text(size=14, face="bold"), 
      plot.title = element_text(size=16), 
      legend.position = "right", 
      legend.title = element_text(size = 14), 
      legend.text = element_text(size = 14)) +
    ggtitle(main) +
    scale_fill_manual(name = "Groups", values = group_colors, 
      breaks = names(group_colors)) +
    xlab("Features") +
    ylab("Proportions")
  
  if(!is.null(prop_fit)){
    ggp <- ggp + 
      geom_point(data = prop_fit, 
        aes_string(x = "feature_id", y = "proportion", 
          group = "sample_id", fill = "group"), 
        position = position_dodge(width = 0.9), size = 3, shape = 23, 
        alpha = 0.75)
  }
  
  return(ggp)
}    


dm_plotProportions_boxplot1 <- function(prop_samp, prop_fit = NULL, 
  main = NULL, group_colors){
  
  ## Plotting
  ggp <- ggplot() +
    geom_jitter(data = prop_samp, aes_string(x = "feature_id", 
      y = "proportion", fill = "group", colour = "group"), 
      position = position_jitterdodge(), 
      alpha = 0.9, size = 2, show.legend = FALSE, na.rm = TRUE) +
    geom_boxplot(data = prop_samp, aes_string(x = "feature_id", 
      y = "proportion", colour = "group", fill = "group"), 
      outlier.size = NA, alpha = 0.4, lwd = 0.5) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
      axis.text=element_text(size=16), 
      axis.title=element_text(size=14, face="bold"), 
      plot.title = element_text(size=16), 
      legend.position = "right", 
      legend.title = element_text(size = 14), 
      legend.text = element_text(size = 14)) +
    ggtitle(main) +     
    scale_fill_manual(name = "Groups", values = group_colors, 
      breaks = names(group_colors)) +
    scale_colour_manual(name = "Groups", values = group_colors, 
      breaks = names(group_colors)) +
    xlab("Features") +
    ylab("Proportions")
  
  if(!is.null(prop_fit)){
    ggp <- ggp + 
      geom_point(data = prop_fit, aes_string(x = "feature_id", 
        y = "proportion", fill = "group"), 
        position = position_jitterdodge(jitter.width = 0), 
        size = 3, shape = 23) +
      guides(colour=FALSE)
  }
  
  return(ggp)
  
}




dm_plotProportions_lineplot <- function(prop_samp, prop_fit = NULL, 
  main = NULL, group_colors){
  
  ## Plotting
  ggp <- ggplot() +
    geom_line(data = prop_samp, aes_string(x = "feature_id", 
      y = "proportion", group = "sample_id",  colour = "group"), 
      size = 1.1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
      axis.text=element_text(size=16), 
      axis.title = element_text(size=14, face="bold"), 
      plot.title = element_text(size=16), 
      legend.position = "right", 
      legend.title = element_text(size = 14), 
      legend.text = element_text(size = 14)) +
    ggtitle(main) +
    scale_fill_manual(name = "Groups", values = group_colors, 
      breaks = names(group_colors)) +
    scale_colour_manual(name = "Groups", values = group_colors, 
      breaks = names(group_colors)) +
    xlab("Features") +
    ylab("Proportions")
  
  if(!is.null(prop_fit)){
    ggp <- ggp + 
      geom_point(data = prop_fit, aes_string(x = "feature_id", 
        y = "proportion", group = "group", fill = "group"), 
        size = 3, shape = 23) +
      guides(colour=FALSE)
  }
  
  return(ggp)
  
}



dm_plotProportions_boxplot2 <- function(prop_samp, prop_fit = NULL, 
  main = NULL, feature_colors){
  
  ## Plotting
  ggp <- ggplot() +
    geom_boxplot(data = prop_samp, aes_string(x = "group", 
      y = "proportion", fill = "feature_id"), 
      outlier.size = NA) + 
    geom_jitter(data = prop_samp, aes_string(x = "group", 
      y = "proportion", fill = "feature_id"),
      position = position_jitterdodge(), 
      shape = 21, show.legend = FALSE, na.rm = TRUE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
      axis.text=element_text(size=14), 
      axis.title=element_text(size=14, face="bold"), 
      plot.title = element_text(size=14), 
      panel.grid.major = element_blank(), 
      legend.title = element_text(size = 14), 
      legend.text = element_text(size = 14)) +
    geom_vline(xintercept = seq(1, nlevels(prop_samp$group) - 1, 1) + 0.5, 
      color = "gray90") +
    ggtitle(main) +     
    scale_fill_manual(name = "Features", values = feature_colors) +
    guides(fill = guide_legend(nrow = 20)) +
    xlab("Groups") +
    ylab("Proportions")
  
  if(!is.null(prop_fit)){
    ggp <- ggp + 
      geom_point(data = prop_fit, 
        aes_string(x = "group", y = "proportion", fill = "feature_id"), 
        position = position_jitterdodge(jitter.width = 0), 
        size = 3, shape = 23, colour = "black")
  }
  
  return(ggp)
  
}


dm_plotProportions_ribbonplot <- function(prop_fit, 
  main = NULL, feature_colors){
  
  prop_fit_order <- prop_fit[order(prop_fit$feature_id, decreasing = TRUE), ]
  prop_fit_order <- prop_fit_order[order(prop_fit_order$group), ]
  breaks <- unique(prop_fit_order$feature_id)
  width  <- 0.5
  
  ## Get ribbons
  gr <- list()
  
  for (i in 1:(nlevels(prop_fit$group) - 1)){
    # i = 1
    
    prop_fit_ribbon <- 
      prop_fit_order[prop_fit_order$group %in% 
          levels(prop_fit_order$group )[c(i, i+1)], ]
    
    prop_fit_ribbon$group <- factor(prop_fit_ribbon$group)
    prop_fit_ribbon$cumsum <- 
      matrix(t(aggregate(prop_fit_ribbon[,"proportion"], 
        by = list(group = prop_fit_ribbon$group), cumsum)[, -1]), 
        ncol = 1)
    prop_fit_ribbon$offset <- 
      c(width/2, -width/2)[as.numeric(prop_fit_ribbon$group)]
    prop_fit_ribbon$xid <- i - 1
    prop_fit_ribbon$x <- as.numeric(prop_fit_ribbon$group) + 
      prop_fit_ribbon$offset + prop_fit_ribbon$xid
    prop_fit_ribbon$ymin <- prop_fit_ribbon$cumsum - 
      prop_fit_ribbon$proportion
    prop_fit_ribbon$ymax <- prop_fit_ribbon$cumsum
    
    gr[[i]] <- geom_ribbon(data = prop_fit_ribbon, 
      aes_string(x = "x", ymin = "ymin", ymax = "ymax", 
        group = "feature_id", fill = "feature_id"), 
      alpha = 0.3) 
    
  }
  
  ## Plotting
  ggp <- ggplot() +
    geom_bar(data = prop_fit_order, 
      aes_string(x = "group", y = "proportion", fill = "feature_id"), 
      stat = "identity", width = width, position="stack") + gr +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
      axis.text=element_text(size=16), 
      axis.title=element_text(size=14, face="bold"), 
      plot.title = element_text(size=16), 
      legend.title = element_text(size = 14), 
      legend.text = element_text(size = 14)) +
    ggtitle(main) +    
    coord_cartesian(ylim = c(-0.1, 1.1)) + 
    coord_cartesian(ylim = c(-0.1, 1.1)) +
    scale_fill_manual(name = "Features", values = feature_colors, 
      breaks = breaks) +
    guides(fill = guide_legend(nrow = 20)) +
    xlab("Groups") +
    ylab("Estimated proportions") 
  
  return(ggp)
  
}

#' Plot feature proportions
#' 
#' Plot observed and/or estimated feature proportions.
#' 
#' @param counts Matrix with rows corresponding to features and columns 
#'   corresponding to samples. Row names are used as labels on the plot.
#' @param group Factor that groups samples into conditions.
#' @param md Data frame with additional sample information.
#' @param fit_full Matrix of estimated proportions with rows corresponding to 
#'   features and columns corresponding to samples. If \code{NULL}, nothing is
#'   plotted.
#' @param main Character vector with main title for the plot. If \code{NULL}, 
#'   nothing is plotted.
#' @param plot_type Character defining the type of the plot produced. Possible 
#'   values \code{"barplot"}, \code{"boxplot1"}, \code{"boxplot2"}, 
#'   \code{"lineplot"}, \code{"ribbonplot"}.
#' @param order_features Logical. Whether to plot the features ordered by their 
#'   expression.
#' @param order_samples Logical. Whether to plot the samples ordered by the 
#'   group variable. If \code{FALSE} order from the \code{sample(x)} is kept.
#' @param group_colors Character vector with colors for each group.
#' @param feature_colors Character vector with colors for each feature.
#'   
#' @return \code{ggplot} object with the observed and/or estimated with 
#'   Dirichlet-multinomial model feature ratios. Estimated proportions are 
#'   marked with diamond shapes.
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes_string theme_bw xlab ylab theme element_text 
#'   coord_cartesian geom_text ggtitle geom_bar scale_fill_manual geom_point 
#'   geom_jitter position_dodge position_jitterdodge geom_boxplot 
#'   scale_colour_manual scale_colour_manual guides element_blank geom_vline 
#'   scale_x_discrete guide_legend geom_line geom_ribbon
#' @importFrom stats aggregate median

dm_plotProportions <- function(counts, group, md = NULL, fit_full = NULL, 
  main = NULL, plot_type = "boxplot1", 
  order_features = TRUE, order_samples = TRUE,
  group_colors = NULL, feature_colors = NULL){
  
  ## Calculate observed proportions
  proportions <- prop.table(counts, 2)
  proportions[proportions == "NaN"] <- NA
  
  prop_samp <- data.frame(feature_id = rownames(proportions), proportions, 
    stringsAsFactors = FALSE) 
  
  prop_fit <- NULL
  
  if(!is.null(fit_full))
    prop_fit <- data.frame(feature_id = rownames(fit_full), fit_full, 
      stringsAsFactors = FALSE)
  
  ## Order transcipts by decreasing proportion 
  if(order_features){
    oo <- order(apply(aggregate(t(prop_samp[, -1]), 
      by = list(group = group), median)[, -1], 2, max), decreasing = TRUE)
    feature_levels <- rownames(prop_samp)[oo]  
  }else{
    feature_levels <- rownames(counts)
  }
  
  ## Order samples by group
  if(order_samples){
    o <- order(group)
    sample_levels <- colnames(counts)[o]
  }else{
    sample_levels <- colnames(counts)
  }
  
  ## Melt prop_samp
  prop_samp <- melt(prop_samp, id.vars = "feature_id", 
    variable.name = "sample_id", value.name = "proportion", 
    factorsAsStrings = FALSE)
  
  prop_samp$feature_id <- factor(prop_samp$feature_id, levels = feature_levels)
  prop_samp$group <- rep(group, each = nrow(counts))
  prop_samp$sample_id <- factor(prop_samp$sample_id, levels = sample_levels)
  
  ## Add extra info from md about samples
  if(!is.null(md)){
    mm <- match(prop_samp$sample_id, md$sample_id)
    for(i in setdiff(colnames(md), c("sample_id", "group"))){
      prop_samp[, i] <- md[mm, i]
    }
  }
  
  ## Melt prop_fit
  if(!is.null(prop_fit)){
    
    prop_fit <- melt(prop_fit, id.vars = "feature_id", 
      variable.name = "sample_id", value.name = "proportion", 
      factorsAsStrings = FALSE)
    prop_fit$feature_id <- factor(prop_fit$feature_id, levels = feature_levels)
    prop_fit$group <- rep(group, each = nrow(fit_full))
    prop_fit$sample_id <- factor(prop_fit$sample_id, levels = sample_levels)
    
    ## Add extra info from md about samples
    if(!is.null(md)){
      mm <- match(prop_fit$sample_id, md$sample_id)
      for(i in setdiff(colnames(md), c("sample_id", "group"))){
        prop_fit[, i] <- md[mm, i]
      }
    }
    
  }
  
  ## Prepare colors for groups
  if(is.null(group_colors)){
    group_colors <- colorb(nlevels(group))
  }
  names(group_colors) <- levels(group)
  
  ## Prepare colors for features
  if(is.null(feature_colors)){
    feature_colors <- colorb(nrow(counts))
  }
  names(feature_colors) <- rownames(counts)
  
  
  switch(plot_type, 
    
    barplot = {
      
      dm_plotProportions_barplot(prop_samp = prop_samp, prop_fit = prop_fit, 
        main = main, group_colors = group_colors)
    },
    
    lineplot = {
      
      dm_plotProportions_lineplot(prop_samp = prop_samp, prop_fit = prop_fit, 
        main = main, group_colors = group_colors)
      
    },
    
    boxplot1 = {
      
      dm_plotProportions_boxplot1(prop_samp = prop_samp, prop_fit = prop_fit, 
        main = main, group_colors = group_colors)
      
    },
    
    boxplot2 = {
      
      dm_plotProportions_boxplot2(prop_samp = prop_samp, prop_fit = prop_fit, 
        main = main, feature_colors = feature_colors)
      
    },
    
    ribbonplot = {
      
      if(!is.null(prop_fit)){
        
        keep <- !duplicated(prop_fit[, c("feature_id", "proportion", "group")])
        
        prop_fit <- prop_fit[keep, , drop = FALSE]
        
        if(nlevels(factor(prop_fit$sample_id)) != 
            nlevels(factor(prop_fit$group))){
          message("Ribbonplot can not be generated.")
        }else{
          dm_plotProportions_ribbonplot(prop_fit,
            main = main, feature_colors = feature_colors)
        }
        
      }else{
        message("Ribbonplot can not be generated.")
      }
      
    }
    
  )
  
}












