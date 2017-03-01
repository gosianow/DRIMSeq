

dmSQTL_CRadjustmentManyGroups_gene <- function(g, counts, genotypes,
  disp, fit, verbose){  
  
  if(verbose >= 2) message(" Gene:", g)
  
  y <- counts[[g]]
  x <- genotypes[[g]]
  
  a <- numeric(nrow(x))
  
  for(i in 1:nrow(x)){
    # i = 2
    
    NAs <- is.na(x[i, ]) | is.na(y[1, ])            
    yy <- y[, !NAs, drop = FALSE]             
    xx <- x[i, !NAs]
    ff <- fit[[g]][[i]][, !NAs, drop = FALSE]
    
    groups <- factor(xx)
    ngroups <- nlevels(groups)
    lgroups <- levels(groups)
    igroups <- lapply(lgroups, function(gr){which(groups == gr)})
    names(igroups) <- lgroups
    # Get the column number of a first occurance of a group level
    figroups <- unlist(lapply(igroups, function(x){x[1]}))
    
    a[i] <- dm_CRadjustmentManyGroups(y = yy, 
      ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
      disp = disp[[g]][i], prop = ff[, figroups, drop = FALSE])
    
  }
  
  # a vector of length #snps for gene g
  return(a)
  
}

#' @importFrom stats model.matrix

dmSQTL_CRadjustmentRegression_gene <- function(g, counts, genotypes,
  group_formula = ~ group, disp, fit, verbose){  
  
  if(verbose >= 2) message(" Gene:", g)
  
  y <- counts[[g]]
  x <- genotypes[[g]]
  
  a <- numeric(nrow(x))
  
  for(i in 1:nrow(x)){
    # i = 2
    
    NAs <- is.na(x[i, ]) | is.na(y[1, ])            
    yy <- y[, !NAs, drop = FALSE]             
    xx <- x[i, !NAs]
    ff <- fit[[g]][[i]][, !NAs, drop = FALSE]
    
    design <- model.matrix(group_formula, data = data.frame(group = xx))
    
    a[i] <- dm_CRadjustmentRegression(y = yy, x = design, 
      disp = disp[[g]][i], prop = ff)
    
  }
  
  # a vector of length #snps for gene g
  return(a)
  
}



dmSQTL_CRadjustment <- function(counts, fit, genotypes, 
  group_formula = ~ group, dispersion, one_way = TRUE,
  verbose = FALSE, BPPARAM = BiocParallel::SerialParam()){
  
  time_start <- Sys.time()
  if(verbose) message("* Calculating Cox-Reid adjustment.. \n")
  
  inds <-  1:length(counts)
  
  # Prepare dispersion
  if(class(dispersion) == "numeric"){ 
    disp <- relist(rep(dispersion, nrow(genotypes)), genotypes@partitioning)
  } else {
    disp <- dispersion
  }
  
  if(one_way){
    
    adj <- BiocParallel::bplapply(inds, 
      dmSQTL_CRadjustmentManyGroups_gene, counts = counts, 
      genotypes = genotypes, disp = disp, fit = fit, 
      verbose = verbose, BPPARAM = BPPARAM)
    
    names(adj) <- names(counts)
    
  }else{
    
    adj <- BiocParallel::bplapply(inds, 
      dmSQTL_CRadjustmentRegression_gene, counts = counts, 
      genotypes = genotypes, group_formula = ~ group, disp = disp, fit = fit, 
      verbose = verbose, BPPARAM = BPPARAM)
    
    names(adj) <- names(counts)
    
  }
  
  time_end <- Sys.time()
  if(verbose >= 2) message("\n")
  if(verbose) message("Took ", round(time_end - time_start, 4), " seconds.\n")
  
  # adj is a list of length G
  return(adj)
  
}












