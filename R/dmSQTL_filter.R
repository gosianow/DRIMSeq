

dmSQTL_filter_genotypes_per_gene <- function(g, counts_new, genotypes, 
  minor_allele_freq){ 
  # g = 1
  
  counts_gene <- counts_new[[g]]
  genotypes_gene <- genotypes[[g]]
  
  ## NA for samples with non expressed genes and missing genotype
  genotypes_gene[, is.na(counts_gene[1,])] <- NA
  genotypes_gene[genotypes_gene == -1] <- NA
  
  ### Keep genotypes with at least minor_allele_freq number of 
  ### variants per group; in other case replace them with NAs
  genotypes_gene <- apply(genotypes_gene, 1, function(x){
    # x <- genotypes_gene[6,]
    
    tt <- table(x)
    
    if( length(tt)==1 )
      return(NULL)
    if( length(tt)==2 ){
      if(any(tt <= minor_allele_freq))
        return(NULL)
      return(x)
    }else{
      if(sum(tt <= minor_allele_freq) >= 2)
        return(NULL)
      x[x == names(tt[tt <= minor_allele_freq])] <- NA
      return(x)
    }    
  })
  
  if(!is.null(genotypes_gene)){
    if(is.list(genotypes_gene))
      genotypes_gene <- do.call(rbind, genotypes_gene)
    else
      genotypes_gene <- t(genotypes_gene)
  }
  
  return(genotypes_gene)
  
}



dmSQTL_filter <- function(counts, genotypes, blocks, samples, 
  min_samps_gene_expr = 70, min_gene_expr = 20, min_samps_feature_expr = 5, 
  min_feature_expr = 20, min_samps_feature_prop = 5, min_feature_prop = 0.05, 
  max_features = Inf, minor_allele_freq = 5, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  ########################################################
  # filtering on counts, put NA for samples with low gene expression
  ########################################################
  
  inds <- which(elementNROWS(counts) > 1)
  
  counts_new <- lapply(inds, function(g){
    # g = 1
    
    expr_features <- counts[[g]]
    
    ### no genes with no expression
    if(sum(expr_features, na.rm = TRUE) == 0)
      return(NULL)
    
    ### genes with min expression
    if(! sum(colSums(expr_features) >= min_gene_expr, na.rm = TRUE) >= 
        min_samps_gene_expr )
      return(NULL)
    
    ### no features with no expression
    features2keep <- rowSums(expr_features > 0, na.rm = TRUE) > 0
    
    ### no genes with one feature
    if(sum(features2keep) <= 1)
      return(NULL)
    
    expr_features <- expr_features[features2keep, , drop = FALSE]
    
    ### features with min expression
    features2keep <- rowSums(expr_features >= min_feature_expr, na.rm = TRUE) >= 
      min_samps_feature_expr
    
    ### no genes with one feature
    if(sum(features2keep) <= 1)
      return(NULL)
    
    expr_features <- expr_features[features2keep, , drop = FALSE]
    
    ### genes with zero expression
    samps2keep <- colSums(expr_features) > 0 & !is.na(expr_features[1, ])
    
    if(sum(samps2keep) < max(1, min_samps_feature_prop))
      return(NULL)
    
    ### consider only samples that have min gene expression, to other assign NAs
    samps2keep <- colSums(expr_features) >= min_gene_expr & 
      !is.na(expr_features[1, ])
    
    if(sum(samps2keep) < max(1, min_samps_feature_prop))
      return(NULL)
    
    prop <- prop.table(expr_features[, samps2keep, drop = FALSE], 2) 
    features2keep <- rowSums(prop >= min_feature_prop) >= 
      min_samps_feature_prop
    
    ### no genes with one feature
    if(sum(features2keep) <= 1)
      return(NULL)

    expr <- expr_features[features2keep, , drop = FALSE] 
    expr[, !samps2keep] <- NA
    
    return(expr)
    
  })
  
  names(counts_new) <- names(counts)[inds]
  NULLs <- !sapply(counts_new, is.null)
  counts_new <- counts_new[NULLs]
  
  if(length(counts_new) == 0)
    stop("!No genes left after filtering!")
  
  counts_new <- MatrixList(counts_new)
  
  ########################################################
  # filtering on genotypes
  ########################################################
  
  genotypes <- genotypes[inds[NULLs]]
  blocks <- blocks[inds[NULLs], ]
  
  genotypes_new <- BiocParallel::bplapply(1:length(counts_new), 
    dmSQTL_filter_genotypes_per_gene, 
    counts_new = counts_new, genotypes = genotypes, 
    minor_allele_freq = minor_allele_freq, BPPARAM = BPPARAM)
  
  names(genotypes_new) <- names(genotypes)
  NULLs <- !sapply(genotypes_new, is.null)
  genotypes_new <- genotypes_new[NULLs]
  
  if(length(genotypes_new) == 0)
    stop("!No SNPs left after filtering!")
  
  genotypes_new <- MatrixList(genotypes_new)
  counts_new <- counts_new[NULLs]
  
  blocks <- blocks[NULLs, ]
  
  ########################################################
  # filtering on blocks
  ########################################################
  
  inds <- 1:length(genotypes_new)
  
  blocks_new <- MatrixList(lapply(inds, function(b){
    # b = 1
    blocks[[b]][blocks[[b]][, "block_id"] %in% rownames(genotypes_new[[b]]), , 
      drop = FALSE]
  }))
  names(blocks_new) <- names(genotypes_new)
  
  data <- new("dmSQTLdata", counts = counts_new, genotypes = genotypes_new, 
    blocks = blocks_new, samples = samples)
  
  return(data)
  
}









