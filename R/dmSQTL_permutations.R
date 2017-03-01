

dmSQTL_permutations_all_genes <- function(x, pvalues, 
  max_nr_perm_cycles = 10, max_nr_min_nr_sign_pval = 1e3, 
  one_way = TRUE, 
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = 0, BPPARAM = BiocParallel::SerialParam()){
	# x dmSQTLfit object
	# pvalues vector of nominal p-values

  nas <- is.na(pvalues)
  pvalues <- pvalues[!nas]
  pvalues <- factor(pvalues)

  sum_sign_pval <- rep(0, length(pvalues))
  nr_perm_tot <- 0
  nr_perm_cycles <- 0
  min_nr_sign_pval <- 0

  while(nr_perm_cycles < max_nr_perm_cycles && 
      min_nr_sign_pval < max_nr_min_nr_sign_pval){
    
    if(verbose)
      message(paste0("** Running permutation cycle number ", 
        nr_perm_cycles + 1 , ".."))

    ### Permute counts for all genes
    n <- ncol(x@counts)
    permutation <- sample(n, n)
    counts_perm <- x@counts[, permutation]
    
    # Fit the DM full model
    fit_full_perm <- dmSQTL_fit(counts = counts_perm, genotypes = x@genotypes, 
    dispersion = x@genewise_dispersion,
    one_way = one_way, group_formula = ~ group,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol,
    return_fit = FALSE, return_coef = FALSE,
    verbose = verbose, BPPARAM = BPPARAM)
    
    # Prepare null (one group) genotypes
    genotypes_null <- x@genotypes
    genotypes_null@unlistData[!is.na(genotypes_null@unlistData)] <- 1
    
    # Fit the DM null model
    fit_null_perm <- dmSQTL_fit(counts = counts_perm, genotypes = genotypes_null, 
      dispersion = x@genewise_dispersion,
      one_way = one_way, group_formula = ~ 1,
      prop_mode = prop_mode, prop_tol = prop_tol, 
      coef_mode = coef_mode, coef_tol = coef_tol,
      return_fit = FALSE, return_coef = FALSE,
      verbose = verbose, BPPARAM = BPPARAM)

    ## Perform the LR test 
    pval_perm <- lapply(1:length(counts_perm), function(g){
      # g = 1
      
      ## Calculate the degrees of freedom
      df <- (nrow(counts_perm[[g]]) - 1) * 
        (apply(x@genotypes[[g]], 1, function(xx) length(unique(xx))) - 1)
      
      out <- dm_LRT(lik_full = fit_full_perm[["lik"]][[g]], 
      	lik_null = fit_null_perm[["lik"]][[g]], 
        df = df, verbose = FALSE)
      
      return(out[, "pvalue"])
      
    })

    pval_perm <- unlist(pval_perm)

    nr_perm <- length(pval_perm)
    nr_perm_tot <- nr_perm_tot + nr_perm
    nr_perm_cycles <- nr_perm_cycles + 1
    
    ### Count how many pval_permuted is lower than pvalues from the model
    nas_perm <- is.na(pval_perm)
    pval_perm <- pval_perm[!nas_perm]
    pval_perm_cut <- cut(pval_perm, c(-1, levels(pvalues), 2), right=FALSE)
    pval_perm_sum <- table(pval_perm_cut)
    pval_perm_cumsum <- cumsum(pval_perm_sum)[-length(pval_perm_sum)]
    names(pval_perm_cumsum) <- levels(pvalues)
    sum_sign_pval <- sum_sign_pval + pval_perm_cumsum[pvalues]
    
    pval_adj <- (sum_sign_pval + 1) / (nr_perm_tot + 1)
    
    min_nr_sign_pval <- min(sum_sign_pval)
    
  }

  pval_out <- rep(NA, length(pvalues))
  pval_out[!nas] <- pval_adj
  
  # pval_out vector of permutation adjusted p-values
  return(pval_out)
  
}


dmSQTL_permutations_per_gene <- function(x, pvalues,
  max_nr_perm = 1e6, max_nr_sign_pval = 1e2,
  one_way = TRUE, 
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = 0, BPPARAM = BiocParallel::SerialParam()){
  # x dmSQTLfit object
  # pvalues list of length G of vector with nominal p-values
  
  pvalues <- lapply(pvalues, function(x){
    pval_tmp <- x[!is.na(x)]
    pval_tmp <- factor(pval_tmp)
    return(pval_tmp)
  })
  
  results_width <- unlist(lapply(pvalues, length))
  nas <- results_width == 0
  genes2permute <- which(!nas)
  
  sum_sign_pval <- vector("list", length(pvalues))
  sum_sign_pval[!nas] <- split(rep(0, sum(results_width)), 
    factor(rep(1:length(results_width), times = results_width)))
  
  nr_perm_tot <- rep(0, length(x@counts))
  nr_perm_tot[nas] <- NA
  
  min_nr_sign_pval <- rep(0, length(x@counts))
  min_nr_sign_pval[nas] <- NA
  
  n <- ncol(x@counts)
  
  while(length(genes2permute) > 0){
    
    if(verbose)
      message(paste0("** ", length(genes2permute), 
        " genes left for permutation.."))

    ### Permute counts for all genes that need additional permutations
    permutation <- sample(n, n)
    counts_perm <- x@counts[genes2permute, permutation]
    genotypes <- x@genotypes[genes2permute, ]
    dispersion <- x@genewise_dispersion[genes2permute]
    
    # Fit the DM full model
    fit_full_perm <- dmSQTL_fit(counts = counts_perm, genotypes = genotypes, 
    dispersion = dispersion,
    one_way = one_way, group_formula = ~ group,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol,
    return_fit = FALSE, return_coef = FALSE,
    verbose = verbose, BPPARAM = BPPARAM)
    
    # Prepare null (one group) genotypes
    genotypes_null <- genotypes
    genotypes_null@unlistData[!is.na(genotypes_null@unlistData)] <- 1
    
    # Fit the DM null model
    fit_null_perm <- dmSQTL_fit(counts = counts_perm, genotypes = genotypes_null, 
      dispersion = dispersion,
      one_way = one_way, group_formula = ~ 1,
      prop_mode = prop_mode, prop_tol = prop_tol, 
      coef_mode = coef_mode, coef_tol = coef_tol,
      return_fit = FALSE, return_coef = FALSE,
      verbose = verbose, BPPARAM = BPPARAM)

    ## Perform the LR test 
    pval_perm <- lapply(1:length(counts_perm), function(g){
      # g = 1
      
      ## Calculate the degrees of freedom
      df <- (nrow(counts_perm[[g]]) - 1) * 
        (apply(genotypes[[g]], 1, function(xx) length(unique(xx))) - 1)
      
      out <- dm_LRT(lik_full = fit_full_perm[["lik"]][[g]], 
      	lik_null = fit_null_perm[["lik"]][[g]], 
        df = df, verbose = FALSE)
      
      return(out[, "pvalue"])
      
    })

    ### Count how many pval_perm is lower than pvalues from the model
    update_nr_sign_pval <- lapply(1:length(pval_perm), function(i){
        # i = 1
        
        pval_perm_gene <- pval_perm[[i]]
        pval_perm_gene <- pval_perm_gene[!is.na(pval_perm_gene)]
        
        pval_perm_cut <- cut(pval_perm_gene, 
          c(-1, levels(pvalues[[genes2permute[i]]]), 2), right=FALSE)

        pval_perm_sum <- table(pval_perm_cut)
        pval_perm_cumsum <- cumsum(pval_perm_sum)[-length(pval_perm_sum)]
        names(pval_perm_cumsum) <- levels(pvalues[[genes2permute[i]]])

        nr_sign_pval <- pval_perm_cumsum[pvalues[[genes2permute[i]]]]
        names(nr_sign_pval) <- NULL

        return(nr_sign_pval)
        
      })

    ### Update values in sum_sign_pval
    for(i in 1:length(update_nr_sign_pval)){
      sum_sign_pval[[genes2permute[i]]] <- sum_sign_pval[[genes2permute[i]]] + 
        update_nr_sign_pval[[i]]
    }
    
    nr_perm <- unlist(lapply(pval_perm, length))
    nr_perm_tot[genes2permute] <- nr_perm_tot[genes2permute] + nr_perm
    
    min_nr_sign_pval[genes2permute] <- 
      unlist(lapply(sum_sign_pval[genes2permute], min))
    
    ### Update genes2permute
    genes2permute <- which(nr_perm_tot < max_nr_perm & 
        min_nr_sign_pval < max_nr_sign_pval)
    
  }

  ### Calculate permutation adjusted p-values
  pval_adj <- lapply(1:length(pvalues), function(i){
      
      pval_tmp <- pvalues[[i]]
      nas <- is.na(pval_tmp)
      
      if(sum(!nas) == 0)
        return(pval_tmp)
      
      pval_tmp[!nas] <- (sum_sign_pval[[i]] + 1) / (nr_perm_tot[i] + 1)
      
      return(pval_tmp)
      
    })
  
  # pval_out list of length G with vectors of permutation adjusted p-values
  return(pval_out)
  
}






