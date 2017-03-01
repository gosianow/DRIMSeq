

dmSQTL_permutations_all_genes <- function(x, fit_null, results, 
  max_nr_perm_cycles = 10, max_nr_min_nr_sign_pval = 1e3, 
  prop_mode, prop_tol, verbose, BPPARAM){
  # results is a data.frame
  
  fit_full <- x@fit_full
  n <- ncol(x@counts)
  
  nr_perm_tot <- 0
  nr_perm_cycles <- 0
  min_nr_sign_pval <- 0
  
  pval <- results$pvalue
  nas <- is.na(pval)
  pval <- pval[!nas]
  pval <- factor(pval)
  sum_sign_pval <- rep(0, length(pval))
  
  pval_perm_all <- matrix(NA, ncol = max_nr_perm_cycles, nrow = nrow(results))
  
  # ds_genes <- results$adj_pvalue < 0.1
  
  while(nr_perm_cycles < max_nr_perm_cycles && 
      min_nr_sign_pval < max_nr_min_nr_sign_pval){
    
    if(verbose)
      message(paste0("** Running cycle number ", nr_perm_cycles + 1 , ".."))
    
    permutation <- sample(n, n)
    
    ### Permute counts for all genes
    counts <- x@counts[, permutation, drop = FALSE]
    
    fit_full_perm <- dmSQTL_fitOneModel(counts = counts, 
      genotypes = x@genotypes, dispersion = slot(x, x@dispersion), 
      model = "full", prop_mode = prop_mode, prop_tol = prop_tol, 
      verbose = verbose, BPPARAM = BPPARAM)
    
    fit_null_perm <- dmSQTL_fitOneModel(counts = counts, 
      genotypes = x@genotypes, dispersion = slot(x, x@dispersion), 
      model = "null", prop_mode = prop_mode, prop_tol = prop_tol, 
      verbose = verbose, BPPARAM = BPPARAM)
    
    results_perm <- dmSQTL_test(fit_full = fit_full_perm, 
      fit_null = fit_null_perm, verbose = verbose, BPPARAM = BPPARAM)
    
    nr_perm <- nrow(results_perm)
    nr_perm_tot <- nr_perm_tot + nr_perm
    nr_perm_cycles <- nr_perm_cycles + 1
    
    
    ### Count how many pval_permuted is lower than pval from the model
    pval_perm_all[, nr_perm_cycles] <- results_perm$pvalue
    nas_perm <- is.na(results_perm$pvalue)
    pval_perm <- results_perm$pvalue[!nas_perm]
    pval_perm_cut <- cut(pval_perm, c(-1, levels(pval), 2), right=FALSE)
    pval_perm_sum <- table(pval_perm_cut)
    pval_perm_cumsum <- cumsum(pval_perm_sum)[-length(pval_perm_sum)]
    names(pval_perm_cumsum) <- levels(pval)
    sum_sign_pval <- sum_sign_pval + pval_perm_cumsum[pval]
    
    pval_adj <- (sum_sign_pval + 1) / (nr_perm_tot + 1)
    
    min_nr_sign_pval <- min(sum_sign_pval)
    
    
  }
  
  
  pval_out <- rep(NA, nrow(results))
  pval_out[!nas] <- pval_adj
  
  return(pval_out)
  
}



dmSQTL_permutations_per_gene <- function(x, fit_null, results, 
  max_nr_perm = 1e6, max_nr_sign_pval = 1e2, prop_mode, prop_tol, 
  verbose = TRUE, BPPARAM){
  # results is a list of data.frames
  
  fit_full <- x@fit_full
  n <- ncol(x@counts)
  
  pval <- lapply(results, function(x){
    pval_tmp <- x$pvalue
    pval_tmp <- pval_tmp[!is.na(pval_tmp)]
    pval_tmp <- factor(pval_tmp)
    return(pval_tmp)
  })
  
  results_width <- unlist(lapply(pval, length))
  nas <- results_width == 0
  genes2permute <- which(!nas)
  
  sum_sign_pval <- vector("list", length(pval))
  sum_sign_pval[!nas] <- split(rep(0, sum(results_width)), 
    factor(rep(1:length(results_width), times = results_width)))
  
  nr_perm_tot <- rep(0, length(x@counts))
  nr_perm_tot[nas] <- NA
  
  min_nr_sign_pval <- rep(0, length(x@counts))
  min_nr_sign_pval[nas] <- NA
  
  
  while(length(genes2permute) > 0){
    
    if(verbose)
      message(paste0("** ", length(genes2permute), 
        " genes left for permutation.."))
    
    permutation <- sample(n, n)
    
    ### Permute counts for all genes that need additional permutations
    counts <- x@counts[genes2permute, permutation]
    genotypes <- x@genotypes[genes2permute, ]
    dispersion <- slot(x, x@dispersion)
    if(is.list(dispersion))
      dispersion <- dispersion[genes2permute]
    
    fit_full_perm <- dmSQTL_fitOneModel(counts = counts, genotypes = genotypes, 
      dispersion = dispersion, model = "full", prop_mode = prop_mode, 
      prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
    fit_null_perm <- dmSQTL_fitOneModel(counts = counts, genotypes = genotypes, 
      dispersion = dispersion, model = "null", prop_mode = prop_mode, 
      prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
    results_perm <- dmSQTL_test(fit_full = fit_full_perm, 
      fit_null = fit_null_perm, return_list = TRUE, 
      verbose = verbose, BPPARAM = BPPARAM)
    
    
    ### Count how many pval_permuted is lower than pval from the model
    update_nr_sign_pval <- lapply(1:length(results_perm), 
      function(i, results_perm, pval, genes2permute){
        # i = 1
        
        pval_perm <- results_perm[[i]]$pvalue
        pval_perm <- pval_perm[!is.na(pval_perm)]
        
        pval_perm_cut <- cut(pval_perm, 
          c(-1, levels(pval[[genes2permute[i]]]), 2), right=FALSE)
        pval_perm_sum <- table(pval_perm_cut)
        
        pval_perm_cumsum <- cumsum(pval_perm_sum)[-length(pval_perm_sum)]
        names(pval_perm_cumsum) <- levels(pval[[genes2permute[i]]])
        nr_sign_pval <- pval_perm_cumsum[pval[[genes2permute[i]]]]
        
        return(nr_sign_pval)
        
      }, results_perm = results_perm, pval = pval, 
      genes2permute = genes2permute)
    
    
    ### Update values in sum_sign_pval
    for(i in 1:length(update_nr_sign_pval)){
      sum_sign_pval[[genes2permute[i]]] <- sum_sign_pval[[genes2permute[i]]] + 
        update_nr_sign_pval[[i]]
    }
    
    
    nr_perm <- unlist(lapply(results_perm, nrow))
    nr_perm_tot[genes2permute] <- nr_perm_tot[genes2permute] + nr_perm
    
    min_nr_sign_pval[genes2permute] <- 
      unlist(lapply(sum_sign_pval[genes2permute], min))
    
    ### Update genes2permute
    genes2permute <- which(nr_perm_tot < max_nr_perm & 
        min_nr_sign_pval < max_nr_sign_pval)
    
  }
  
  
  ### Calculate permutation adjusted p-values
  pval_adj <- lapply(1:length(results), 
    function(i, results, sum_sign_pval, nr_perm_tot){
    
    pval_tmp <- results[[i]]$pvalue
    nas <- is.na(pval_tmp)
    
    if(sum(!nas) == 0)
      return(pval_tmp)
    
    pval_tmp[!nas] <- (sum_sign_pval[[i]] + 1) / (nr_perm_tot[i] + 1)
    
    return(pval_tmp)
    
  }, results = results, sum_sign_pval = sum_sign_pval, 
    nr_perm_tot = nr_perm_tot)
  
  
  pval_out <- unlist(pval_adj)
  
  return(pval_out)
  
}






