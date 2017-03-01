#' @include class_dmSQTLfit.R class_dmDStest.R
NULL

###############################################################################
### dmSQTLtest class
###############################################################################

#' dmSQTLtest object
#' 
#' dmSQTLtest extends the \code{\linkS4class{dmSQTLfit}} class by adding the 
#' null model Dirichlet-multinomial feature proportion estimates and the results
#' of testing for sQTLs. Proportions are calculated for each gene-block pair 
#' from pooled (no grouping into conditions) counts. Result of 
#' \code{\link{dmTest}}.
#' 
#' @return
#' 
#' \itemize{ \item \code{results(x)}: Get a data frame with results. See Slots. 
#' }
#' 
#' @param x dmSQTLtest object.
#' @param ... Other parameters that can be defined by methods using this 
#'   generic.
#'   
#' @slot lik_null List of numeric vectors with the per gene-snp DM null model
#'   likelihoods.
#' @slot results_gene Data frame with the gene-level results including: 
#'   \code{gene_id} - gene IDs, \code{block_id} - block IDs, \code{snp_id} - SNP
#'   IDs, \code{lr} - likelihood ratio statistics based on the DM model, 
#'   \code{df} - degrees of freedom, \code{pvalue} - p-values estimated based on
#'   permutations and \code{adj_pvalue} - Benjamini & Hochberg adjusted 
#'   p-values.
#'   
#' @examples 
#' 
#' #############################
#' ### sQTL analysis
#' #############################
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmSQTLdata}}, \code{\linkS4class{dmSQTLdata}}, 
#'   \code{\linkS4class{dmSQTLdispersion}}, \code{\linkS4class{dmSQTLfit}}
setClass("dmSQTLtest", 
  contains = "dmSQTLfit",
  representation(lik_null = "list",
    results_gene = "data.frame"))


#####################################


setValidity("dmSQTLtest", function(object){
  # has to return TRUE when valid object!
  
  # TODO: Add more checks
  
  if(!length(object@counts) == length(object@lik_null))
    return("Different number of genes in 'counts' and 'lik_null'")
  
  return(TRUE)
  
})


###############################################################################
### show and accessing methods
###############################################################################


#' @rdname dmSQTLtest-class
#' @export
setMethod("results", "dmSQTLtest", function(x) x@results_gene)



# -----------------------------------------------------------------------------

setMethod("show", "dmSQTLtest", function(object){
  
  callNextMethod(object)
  
  cat("* data accessors: results()\n")
  
})


###############################################################################
### dmTest
###############################################################################

#' @param permutation_mode Character specifying which permutation scheme to
#'   apply for p-value calculation. When equal to \code{"all_genes"}, null 
#'   distribution of p-values is calculated from all genes and the maximum 
#'   number of permutation cycles is 10. When  \code{permutation_mode =
#'   "per_gene"}, null distribution of p-values is calculated for each gene
#'   separately based on permutations of this individual gene. The latter
#'   approach may take a lot of computational time. We suggest using the first
#'   option.
#' @rdname dmTest
#' @export
#' @importFrom utils relist
setMethod("dmTest", "dmSQTLfit", function(x, 
  permutation_mode = "all_genes", one_way = TRUE, 
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = 0, BPPARAM = BiocParallel::SerialParam()){
  
  # Check parameters
  stopifnot(permutation_mode %in% c("all_genes", "per_gene"))
  stopifnot(is.logical(one_way))
  
  stopifnot(length(prop_mode) == 1)
  stopifnot(prop_mode %in% c("constrOptim"))
  stopifnot(length(prop_tol) == 1)
  stopifnot(is.numeric(prop_tol) && prop_tol > 0)
  
  stopifnot(length(coef_mode) == 1)
  stopifnot(coef_mode %in% c("optim", "nlminb", "Rcgmin"))
  stopifnot(length(coef_tol) == 1)
  stopifnot(is.numeric(coef_tol) && coef_tol > 0)
  
  stopifnot(verbose %in% 0:2)
  
  # Prepare null (one group) genotypes
  genotypes_null <- x@genotypes
  genotypes_null@unlistData[!is.na(genotypes_null@unlistData)] <- 1
  
  # Fit the DM null model
  fit0 <- dmSQTL_fit(counts = x@counts, genotypes = genotypes_null, 
    dispersion = x@genewise_dispersion,
    one_way = one_way, group_formula = ~ 1,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol,
    return_fit = FALSE, return_coef = FALSE,
    verbose = verbose, BPPARAM = BPPARAM)
  
  ## Perform the LR test 
  results_list <- lapply(1:length(x@counts), function(g){
    # g = 1
    
    ## Calculate the degrees of freedom
    df <- (nrow(x@counts[[g]]) - 1) * 
      (apply(x@genotypes[[g]], 1, function(x) length(unique(x))) - 1)
    
    out <- dm_LRT(lik_full = x@lik_full[[g]], lik_null = fit0[["lik"]][[g]], 
      df = df, verbose = FALSE)
    
    return(out)
    
  })
  
  if(verbose)
    message("\n** Running permutations..\n")
  
  ### Calculate adjusted p-values using permutations
  switch(permutation_mode, 
    all_genes = {
      ## P-value for a gene computed using all the permutations
      
      pvalues <- unlist(lapply(results_list, function(x) x[, "pvalue"]))
      
      pval_adj_perm <- dmSQTL_permutations_all_genes(x = x, pvalues = pvalues, 
        max_nr_perm_cycles = 10, max_nr_min_nr_sign_pval = 1e3, 
        one_way = one_way,
        prop_mode = prop_mode, prop_tol = prop_tol, 
        coef_mode = coef_mode, coef_tol = coef_tol,
        verbose = verbose, BPPARAM = BPPARAM)
      
      pval_adj_perm <- relist(pval_adj_perm, x@lik_full)
      
    },
    
    per_gene = {
      ## P-value for a gene computed using permutations of that gene 
      
      pvalues <- lapply(results_list, function(x) x[, "pvalue"])
      
      pval_adj_perm <- dmSQTL_permutations_per_gene(x = x, pvalues = pvalues, 
        max_nr_perm = 1e6, max_nr_sign_pval = 1e2, 
        one_way = one_way,
        prop_mode = prop_mode, prop_tol = prop_tol, 
        coef_mode = coef_mode, coef_tol = coef_tol,
        verbose = verbose, BPPARAM = BPPARAM)
      
    }
    
  )
  
  pval_adj_perm_BH <- relist(p.adjust(unlist(pval_adj_perm), method="BH"), 
    pval_adj_perm) 
  
  inds <- 1:length(results_list)
  
  for(i in inds){
    results_list[[i]][, "pvalue"] <- pval_adj_perm[[i]]
    results_list[[i]][, "adj_pvalue"] <- pval_adj_perm_BH[[i]]
  }
  
  gene_ids <- names(x@blocks)
  
  ## Output the original SNPs
  results_new <- lapply(inds, function(i){
    # i = 1
    
    mm <- match(x@blocks[[i]][, "block_id"], rownames(x@genotypes[[i]]))
    
    out <- data.frame(gene_id = gene_ids[i], x@blocks[[i]], 
      results_list[[i]][mm, ], stringsAsFactors = FALSE)
    
    return(out)
    
  })
  
  results_new <- do.call(rbind, results_new)
  
  return(new("dmSQTLtest", lik_null = fit0[["lik"]], results_gene = results_new,
    lik_full = x@lik_full, fit_full = x@fit_full,
    mean_expression = x@mean_expression, 
    common_dispersion = x@common_dispersion, 
    genewise_dispersion = x@genewise_dispersion, 
    counts = x@counts, genotypes = x@genotypes,
    blocks = x@blocks, samples = x@samples))
  
})


###############################################################################
### plotPValues
###############################################################################

#' @rdname plotPValues
#' @export
setMethod("plotPValues", "dmSQTLtest", function(x){
  
  ### Plot p-values for unique blocks (not SNPs)
  
  keep <- !duplicated(x@results_gene[, c("gene_id", "block_id"), drop = FALSE])
  
  ggp <- dm_plotPValues(pvalues = x@results_gene[keep, "pvalue"])

  return(ggp)  
  
})

























