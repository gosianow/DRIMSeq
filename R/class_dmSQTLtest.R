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

#' @param permutations Character specifying which permutation scheme to apply
#'   for p-value calculation. When equal to \code{"all_genes"}, null
#'   distribution of p-values is calculated from all genes and the maximum
#'   number of permutation cycles is 10. When  \code{permutations = "per_gene"},
#'   null distribution of p-values is calculated for each gene separately based
#'   on permutations of this individual gene. The latter approach may take a lot
#'   of computational time. We suggest using the first option.
#' @rdname dmTest
#' @export
setMethod("dmTest", "dmSQTLfit", function(x, 
  permutations = "all_genes", one_way = TRUE, 
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = 0, BPPARAM = BiocParallel::SerialParam()){
  
  # Check parameters
  stopifnot(permutations %in% c("all_genes", "per_gene"))
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
  
  results_list <- lapply(1:length(x@counts), function(g){
    # g = 1
    
    df <- (nrow(x@counts[[g]]) - 1) * 
      (apply(x@genotypes[[g]], 1, function(x) length(unique(x))) - 1)
    
    out <- dm_LRT(lik_full = x@lik_full[[g]], lik_null = fit0[["lik"]][[g]], 
      df = df, verbose = FALSE)
    
    return(out)
    
  })
  
  
  if(verbose)
    message("\n** Running permutations..\n")
  
  ### TO DO:
  
  ### Calculate adjusted p-values using permutations
  switch(permutations, 
    all_genes = {
      # P-value for a gene computed using all the permutations
      
      pval_adj_perm <- dmSQTL_permutations_all_genes(x = x, fit_null = fit_null, 
        results = results, max_nr_perm_cycles = 10, max_nr_min_nr_sign_pval = 1e3, 
        prop_mode = prop_mode, prop_tol = prop_tol, 
        verbose = verbose, BPPARAM = BPPARAM)
      
    },
    per_gene = {
      # P-value for a gene computed using permutations of that gene 
      
      pval_adj_perm <- dmSQTL_permutations_per_gene(x = x, fit_null = fit_null, 
        results = results, max_nr_perm = 1e6, max_nr_sign_pval = 1e2, 
        prop_mode = prop_mode, prop_tol = prop_tol, 
        verbose = verbose, BPPARAM = BPPARAM)
      
    }
  )
  
  
  results$pvalue <- pval_adj_perm
  results$adj_pvalue <- p.adjust(pval_adj_perm, method="BH")
  
  # Output the original SNPs
  colnames(results)[colnames(results) == "snp_id"] <- "block_id" 
  results_spl <- split(results, factor(results$gene_id, 
    levels = names(x@blocks)))
  inds <- 1:length(results_spl)
  
  results_new <- lapply(inds, function(i){
    # i = 1
    
    res <- results_spl[[i]]
    blo <- x@blocks[[i]]
    matching <- match(blo[, "block_id"], res[, "block_id"])
    snp_id <- blo[, "snp_id"]
    res_new <- cbind(res[matching, c("gene_id", "block_id")], snp_id, 
      res[matching, c("lr", "df", "pvalue", "adj_pvalue")])
    
    return(res_new)
    
  })
  
  results_new <- do.call(rbind, results_new)
  
  return(new("dmSQTLtest", fit_null = fit_null, results = results_new, 
    dispersion = x@dispersion, fit_full = x@fit_full, 
    mean_expression = x@mean_expression, 
    common_dispersion = x@common_dispersion, 
    genewise_dispersion = x@genewise_dispersion, counts = x@counts, 
    genotypes = x@genotypes, blocks = x@blocks, samples = x@samples))
  
  
})


###############################################################################
### plotPValues
###############################################################################

#' @rdname plotPValues
#' @export
#' @importFrom grDevices pdf dev.off
setMethod("plotPValues", "dmSQTLtest", function(x, out_dir = NULL){
  
  ### Plot p-values for unique blocks (not SNPs)
  ggp <- dm_plotPvalues(pvalues = unique(x@results[, c("gene_id", "block_id", 
    "pvalue"), drop = FALSE])[, "pvalue"])
  
  if(!is.null(out_dir)){
    pdf(paste0(out_dir, "hist_pvalues.pdf"))
    print(ggp)
    dev.off()
  }else{
    return(ggp)  
  }
  
})





