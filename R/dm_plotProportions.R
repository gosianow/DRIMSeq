

dm_plotProportions_barplot <- function(counts, group, prop_full = NULL, 
  main = NULL, order = TRUE, group_colors = NULL){
  
  labels <- labels_org <- factor(rownames(counts), levels = rownames(counts))
  group_counts <- table(group)
  
  proportions <- prop.table(counts, 2)
  proportions[proportions == "NaN"] <- NA
  
  prop_samp <- data.frame(feature_id = labels, proportions, 
    stringsAsFactors = FALSE) 
  
  if(!is.null(prop_full))
    prop_est_full <- data.frame(feature_id = labels, prop_full, 
      stringsAsFactors = FALSE)
  
  #### order transcipts by decreasing proportion 
  if(order){
    labels <- labels[order(apply(aggregate(t(prop_samp[, -1]), 
      by = list(group = group), median)[, -1], 2, max), decreasing = TRUE)]  
  }
  
  prop_samp <- melt(prop_samp, id.vars = "feature_id", 
    variable.name = "sample_id", value.name = "proportion")  
  prop_samp$feature_id <- factor(prop_samp$feature_id, levels = labels)
  prop_samp$sample_id <- factor(prop_samp$sample_id)
  prop_samp$group <- rep(group, each = length(labels))
  
  if(!is.null(prop_full)){
    prop_est_full <- melt(prop_est_full, id.vars = "feature_id", 
      variable.name = "group", value.name = "proportion")
    prop_est_full$feature_id <- factor(prop_est_full$feature_id, 
      levels = labels)
    prop_est_full$group <- factor(rep(levels(group), each = length(labels)), 
      levels = levels(group))
  }
  
  if(is.null(group_colors))
    values <- colorb(nlevels(group))
  else
    values <- group_colors
  names(values) <- levels(group)
  
  order_prop_samp <- order(prop_samp$group, prop_samp$sample_id)
  prop_samp$sample_id <- factor(prop_samp$sample_id, 
    levels = unique(prop_samp$sample_id[order_prop_samp]))
  width = 0.9
  
  if(!is.null(prop_full)){
    prop_est_full$xid <- as.numeric(prop_est_full$feature_id)
    prop_est_full$group_prop <- prop_est_full$group
    levels(prop_est_full$group_prop) <- group_counts/sum(group_counts)
    prop_est_full$group_prop <- 
      as.numeric(as.character(prop_est_full$group_prop))
    prop_est_full$group_cumsum <- prop_est_full$group
    levels(prop_est_full$group_cumsum) <- 
      cumsum(group_counts)/sum(group_counts)
    prop_est_full$group_cumsum <- 
      as.numeric(as.character(prop_est_full$group_cumsum))
    prop_est_full$x <- prop_est_full$xid - width/2 + 
      prop_est_full$group_cumsum * width - prop_est_full$group_prop/2
  }
  
  ggp <- ggplot() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
      axis.text=element_text(size=16), 
      axis.title = element_text(size=14, face="bold"), 
      plot.title = element_text(size=16), 
      legend.position = "right", 
      legend.title = element_text(size = 14), 
      legend.text = element_text(size = 14)) +
    ggtitle(main) +
    geom_bar(data = prop_samp, 
      aes_string(x = "feature_id", y = "proportion", group = "sample_id", 
        fill = "group"), 
      stat = "identity", position = position_dodge(width = width)) +
    scale_fill_manual(name = "Groups", values = values, 
      breaks = names(values)) +
    xlab("Features") +
    ylab("Proportions")
  
  if(!is.null(prop_full)){
    ggp <- ggp + 
      geom_point(data = prop_est_full, 
        aes_string(x = "x", y = "proportion", fill = "group"), 
        size = 3, shape = 23)
  }
  
  return(ggp)
}    


dm_plotProportions_boxplot1 <- function(counts, group, prop_full = NULL, 
  main = NULL, order = TRUE, group_colors = NULL){
  
  labels <- labels_org <- factor(rownames(counts), levels = rownames(counts))
  group_counts <- table(group)
  
  proportions <- prop.table(counts, 2)
  proportions[proportions == "NaN"] <- NA
  
  prop_samp <- data.frame(feature_id = labels, proportions, 
    stringsAsFactors = FALSE) 
  
  if(!is.null(prop_full))
    prop_est_full <- data.frame(feature_id = labels, prop_full, 
      stringsAsFactors = FALSE)
  
  #### order transcipts by decreasing proportion 
  if(order){
    labels <- labels[order(apply(aggregate(t(prop_samp[, -1]), 
      by = list(group = group), median)[, -1], 2, max), decreasing = TRUE)]  
  }
  
  prop_samp <- melt(prop_samp, id.vars = "feature_id", 
    variable.name = "sample_id", value.name = "proportion")  
  prop_samp$feature_id <- factor(prop_samp$feature_id, levels = labels)
  prop_samp$sample_id <- factor(prop_samp$sample_id)
  prop_samp$group <- rep(group, each = length(labels))
  
  if(!is.null(prop_full)){
    prop_est_full <- melt(prop_est_full, id.vars = "feature_id", 
      variable.name = "group", value.name = "proportion")
    prop_est_full$feature_id <- factor(prop_est_full$feature_id, 
      levels = labels)
    prop_est_full$group <- factor(rep(levels(group), each = length(labels)), 
      levels = levels(group))
  }
  
  ### box plots with points
  if(is.null(group_colors))
    values <- colorb(nlevels(group))
  else
    values <- group_colors
  names(values) <- levels(group)
  
  ggp <- ggplot() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
      axis.text=element_text(size=16), 
      axis.title=element_text(size=14, face="bold"), 
      plot.title = element_text(size=16), 
      legend.position = "right", 
      legend.title = element_text(size = 14), 
      legend.text = element_text(size = 14)) +
    ggtitle(main) +     
    geom_jitter(data = prop_samp, 
      aes_string(x = "feature_id", y = "proportion", fill = "group", 
        colour = "group"), 
      position = position_jitterdodge(dodge.width = 0.75), 
      alpha = 0.5, size = 2, show.legend = FALSE, na.rm = TRUE) +
    geom_boxplot(data = prop_samp, 
      aes_string(x = "feature_id", y = "proportion", colour = "group", 
        fill = "group"), 
      outlier.size = NA, alpha = 0.2, lwd = 0.5) +
    scale_fill_manual(name = "Groups", values = values, 
      breaks = names(values)) +
    scale_colour_manual(name = "Groups", values = values, 
      breaks = names(values)) +
    xlab("Features") +
    ylab("Proportions")
  
  if(!is.null(prop_full)){
    ggp <- ggp + 
      geom_point(data = prop_est_full, 
        aes_string(x = "feature_id", y = "proportion", fill = "group"), 
        position = position_jitterdodge(jitter.width = 0, 
          jitter.height = 0), 
        size = 3, shape = 23) +
      guides(colour=FALSE)
  }
  
  return(ggp)
  
}




dm_plotProportions_lineplot <- function(counts, group, prop_full = NULL, 
  main = NULL, order = TRUE, group_colors = NULL){
  
  labels <- labels_org <- factor(rownames(counts), levels = rownames(counts))
  group_counts <- table(group)
  
  proportions <- prop.table(counts, 2)
  proportions[proportions == "NaN"] <- NA
  
  prop_samp <- data.frame(feature_id = labels, proportions, 
    stringsAsFactors = FALSE) 
  
  if(!is.null(prop_full))
    prop_est_full <- data.frame(feature_id = labels, prop_full, 
      stringsAsFactors = FALSE)
  
  #### order transcipts by decreasing proportion 
  if(order){
    labels <- labels[order(apply(aggregate(t(prop_samp[, -1]), 
      by = list(group = group), median)[, -1], 2, max), decreasing = TRUE)]  
  }
  
  prop_samp <- melt(prop_samp, id.vars = "feature_id", 
    variable.name = "sample_id", value.name = "proportion")  
  prop_samp$feature_id <- factor(prop_samp$feature_id, levels = labels)
  prop_samp$sample_id <- factor(prop_samp$sample_id)
  prop_samp$group <- rep(group, each = length(labels))
  
  if(!is.null(prop_full)){
    prop_est_full <- melt(prop_est_full, id.vars = "feature_id", 
      variable.name = "group", value.name = "proportion")
    prop_est_full$feature_id <- factor(prop_est_full$feature_id, 
      levels = labels)
    prop_est_full$group <- factor(rep(levels(group), each = length(labels)), 
      levels = levels(group))
  }
  
  ### line plots
  if(is.null(group_colors))
    values <- colorb(nlevels(group))
  else
    values <- group_colors
  names(values) <- levels(group)
  
  ggp <- ggplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
      axis.text=element_text(size=16), 
      axis.title = element_text(size=14, face="bold"), 
      plot.title = element_text(size=16), 
      legend.position = "right", 
      legend.title = element_text(size = 14), 
      legend.text = element_text(size = 14)) +
    ggtitle(main) +
    geom_line(data = prop_samp, 
      aes_string(x = "feature_id", y = "proportion", group = "sample_id", 
        colour = "group"), 
      size = 1.1) +
    scale_fill_manual(name = "Groups", values = values, 
      breaks = names(values)) +
    scale_colour_manual(name = "Groups", values = values, 
      breaks = names(values)) +
    xlab("Features") +
    ylab("Proportions")
  
  if(!is.null(prop_full)){
    ggp <- ggp + 
      geom_point(data = prop_est_full, 
        aes_string(x = "feature_id", y = "proportion", group = "group", 
          fill = "group"), 
        size = 3, shape = 23) +
      guides(colour=FALSE)
  }
  
  return(ggp)
  
}



dm_plotProportions_boxplot2 <- function(counts, group, prop_full = NULL, 
  main = NULL, order = TRUE, feature_colors = NULL){
  
  labels <- labels_org <- factor(rownames(counts), levels = rownames(counts))
  group_counts <- table(group)
  
  proportions <- prop.table(counts, 2)
  proportions[proportions == "NaN"] <- NA
  
  prop_samp <- data.frame(feature_id = labels, proportions, 
    stringsAsFactors = FALSE) 
  
  if(!is.null(prop_full))
    prop_est_full <- data.frame(feature_id = labels, prop_full, 
      stringsAsFactors = FALSE)
  
  #### order transcipts by decreasing proportion 
  if(order){
    labels <- labels[order(apply(aggregate(t(prop_samp[, -1]), 
      by = list(group = group), median)[, -1], 2, max), decreasing = TRUE)]  
  }
  
  prop_samp <- melt(prop_samp, id.vars = "feature_id", 
    variable.name = "sample_id", value.name = "proportion")  
  prop_samp$feature_id <- factor(prop_samp$feature_id, levels = labels)
  prop_samp$sample_id <- factor(prop_samp$sample_id)
  prop_samp$group <- rep(group, each = length(labels))
  
  if(!is.null(prop_full)){
    prop_est_full <- melt(prop_est_full, id.vars = "feature_id", 
      variable.name = "group", value.name = "proportion")
    prop_est_full$feature_id <- factor(prop_est_full$feature_id, 
      levels = labels)
    prop_est_full$group <- factor(rep(levels(group), each = length(labels)), 
      levels = levels(group))
  }
  
  ### box plots per group
  if(is.null(feature_colors))
    values <- colorb(length(labels_org))
  else
    values <- feature_colors
  names(values) <- labels_org
  
  ggp <- ggplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
      axis.text=element_text(size=14), 
      axis.title=element_text(size=14, face="bold"), 
      plot.title = element_text(size=14), 
      panel.grid.major = element_blank(), 
      legend.title = element_text(size = 14), 
      legend.text = element_text(size = 14)) +
    geom_vline(xintercept = seq(1, nlevels(group) - 1, 1) + 0.5, 
      color = "gray90") +
    ggtitle(main) +     
    geom_boxplot(data = prop_samp, 
      aes_string(x = "group", y = "proportion", fill = "feature_id"), 
      width = 1) + 
    scale_fill_manual(name = "Features", values = values) +
    scale_x_discrete(labels = paste0(names(group_counts), " (", 
      group_counts, ")" ), name="") +
    guides(fill = guide_legend(nrow = 20)) +
    xlab("Groups") +
    ylab("Proportions")
  
  if(!is.null(prop_full)){
    ggp <- ggp + 
      geom_point(data = prop_est_full, 
        aes_string(x = "group", y = "proportion", fill = "feature_id"), 
        position = position_jitterdodge(jitter.width = 0, 
          jitter.height = 0, dodge.width = 1), 
        size = 3, shape = 23, colour = "black")
  }
  
  return(ggp)
  
}


dm_plotProportions_ribbonplot <- function(counts, group, prop_full = NULL, 
  main = NULL, order = TRUE, feature_colors = NULL){
  
  labels <- labels_org <- factor(rownames(counts), levels = rownames(counts))
  group_counts <- table(group)
  
  proportions <- prop.table(counts, 2)
  proportions[proportions == "NaN"] <- NA
  
  prop_samp <- data.frame(feature_id = labels, proportions, 
    stringsAsFactors = FALSE) 
  
  if(!is.null(prop_full))
    prop_est_full <- data.frame(feature_id = labels, prop_full, 
      stringsAsFactors = FALSE)
  
  #### order transcipts by decreasing proportion 
  if(order){
    labels <- labels[order(apply(aggregate(t(prop_samp[, -1]), 
      by = list(group = group), median)[, -1], 2, max), decreasing = TRUE)]  
  }
  
  prop_samp <- melt(prop_samp, id.vars = "feature_id", 
    variable.name = "sample_id", value.name = "proportion")  
  prop_samp$feature_id <- factor(prop_samp$feature_id, levels = labels)
  prop_samp$sample_id <- factor(prop_samp$sample_id)
  prop_samp$group <- rep(group, each = length(labels))
  
  if(!is.null(prop_full)){
    prop_est_full <- melt(prop_est_full, id.vars = "feature_id", 
      variable.name = "group", value.name = "proportion")
    prop_est_full$feature_id <- factor(prop_est_full$feature_id, 
      levels = labels)
    prop_est_full$group <- factor(rep(levels(group), each = length(labels)), 
      levels = levels(group))
  }
  
  if(!is.null(prop_full)){
    
    if(is.null(feature_colors))
      values <- colorb(length(labels_org))
    else
      values <- feature_colors
    names(values) <- labels_org
    
    breaks <- labels
    width  <- 0.5
    prop_est_full_order <- prop_est_full[order(prop_est_full$group), ]
    
    if(order == TRUE){
      prop_est_full_order <- prop_est_full[order(prop_est_full$group, 
        prop_est_full$proportion), ]
      breaks = rev(prop_est_full_order[prop_est_full_order$group == 
          levels(prop_est_full_order$group)[1], "feature_id"])
    }
    
    ### get ribbons!!!
    gr <- list()
    
    for (i in 1:(nlevels(group) - 1)){
      # i = 2
      prop_est_full_ribbon <- 
        prop_est_full_order[prop_est_full_order$group %in% 
            levels(prop_est_full_order$group )[c(i, i+1)], ]
      prop_est_full_ribbon$group <- factor(prop_est_full_ribbon$group)
      prop_est_full_ribbon$cumsum <- 
        matrix(t(aggregate(prop_est_full_ribbon[,"proportion"], 
          by = list(group = prop_est_full_ribbon$group), cumsum)[, -1]), 
          ncol = 1)
      prop_est_full_ribbon$offset <- 
        c(width/2, -width/2)[as.numeric(prop_est_full_ribbon$group)]
      prop_est_full_ribbon$xid <- i - 1
      prop_est_full_ribbon$x <- as.numeric(prop_est_full_ribbon$group) + 
        prop_est_full_ribbon$offset + prop_est_full_ribbon$xid
      prop_est_full_ribbon$ymin <- prop_est_full_ribbon$cumsum - 
        prop_est_full_ribbon$proportion
      prop_est_full_ribbon$ymax <- prop_est_full_ribbon$cumsum
      
      gr[[i]] <- geom_ribbon(data = prop_est_full_ribbon, 
        aes_string(x = "x", ymin = "ymin", ymax = "ymax", 
          group = "feature_id", fill = "feature_id"), 
        alpha = 0.3) 
      
    }
    
    ggp <- ggplot() +
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
      scale_fill_manual(name = "Features", values = values, 
        breaks = breaks) +
      scale_x_discrete(labels = paste0(names(group_counts), " (", 
        group_counts, ")" ), name="") +
      guides(fill = guide_legend(nrow = 25)) +
      xlab("Groups") +
      ylab("Estimated proportions") +
      geom_bar(data = prop_est_full_order, 
        aes_string(x = "group", y = "proportion", fill = "feature_id"), 
        stat = "identity", width = width, position="stack") + gr
    
  }
  
  return(ggp)
  
}

#' Plot feature proportions
#' 
#' Plot observed and/or estimated feature proportions.
#' 
#' @param counts Matrix with rows corresponding to features and columns 
#'   corresponding to samples. Row names are used as labels on the plot.
#' @param group Factor that groups samples into conditions.
#' @param prop_full Matrix of estimated proportions with rows corresponding to 
#'   features and columns corresponding to conditions defined by factor 
#'   \code{group}. If \code{NULL}, nothing is plotted.
#' @param main Character vector with main title for the plot. If \code{NULL}, 
#'   nothing is plotted.
#' @param plot_type Character defining the type of the plot produced. Possible 
#'   values \code{"barplot"}, \code{"boxplot1"}, \code{"boxplot2"}, 
#'   \code{"lineplot"}, \code{"ribbonplot"}.
#' @param order Logical. Whether to plot the features ordered by their 
#'   expression.
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

dm_plotProportions <- function(counts, group, prop_full = NULL, 
  main = NULL, plot_type = "boxplot1", order = TRUE,
  group_colors = NULL, feature_colors = NULL){
  
  switch(plot_type, 
    
    barplot = {
      
      dm_plotProportions_barplot(counts = counts, group = group, 
        prop_full = prop_full, main = main, 
        order = order, group_colors = group_colors)
    },
    
    boxplot1 = {
      
      dm_plotProportions_boxplot1(counts = counts, group = group, 
        prop_full = prop_full, main = main, 
        order = order, group_colors = group_colors)
      
    },
    
    boxplot2 = {
      
      dm_plotProportions_boxplot2(counts = counts, group = group, 
        prop_full = prop_full, main = main, 
        order = order, feature_colors = feature_colors)
      
    },
    
    lineplot = {
      
      dm_plotProportions_lineplot(counts = counts, group = group, 
        prop_full = prop_full, main = main, 
        order = order, group_colors = group_colors)
      
    },
    
    ribbonplot = {
      
      dm_plotProportions_ribbonplot(counts = counts, group = group, 
        prop_full = prop_full, main = main, 
        order = order, feature_colors = feature_colors)
    }
    
  )
  
}












