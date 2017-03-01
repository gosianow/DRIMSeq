
dmDS_filter <- function(counts, min_samps_gene_expr = 6, 
  min_gene_expr = 10, min_samps_feature_expr = 3, min_feature_expr = 10, 
  min_samps_feature_prop = 3, min_feature_prop = 0.01){
  
  inds <- which(elementNROWS(counts) > 1)
  
  counts_new <- lapply(inds, function(g){
    # g = 117
    # print(g)
    
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
    
    prop <- prop.table(expr_features[, samps2keep, drop = FALSE], 2) 
    # prop.table(matrix(c(1,0), 2, 1), 2)
    # prop.table(matrix(c(0,0), 2, 1), 2)
    # prop.table(matrix(c(0,0, 1, 0), 2, 2), 2)
    
    ### features with min proportion
    features2keep <- rowSums(prop >= min_feature_prop) >= min_samps_feature_prop
    
    ### no genes with one feature
    if(sum(features2keep) <= 1)
      return(NULL)
    
    expr <- expr_features[features2keep, , drop = FALSE] 
    
    return(expr)
    
  })
  
  names(counts_new) <- names(counts)[inds]
  counts_new <- counts_new[!sapply(counts_new, is.null)]
  
  if(length(counts_new) == 0)
    stop("!No genes left after filtering!")
  
  counts_new <- MatrixList(counts_new)
  
  return(counts_new)
  
}













