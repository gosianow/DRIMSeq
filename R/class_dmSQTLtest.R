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
#' @slot fit_null List of \code{\linkS4class{MatrixList}}. Each of them contains
#'   null proportions, likelihoods and degrees of freedom for all the blocks 
#'   (unique SNPs) assigned to a given gene.
#' @slot results Data frame with \code{gene_id} - gene IDs, \code{block_id} - 
#'   block IDs, \code{snp_id} - SNP IDs, \code{lr} - likelihood ratio 
#'   statistics, \code{df} - degrees of freedom, \code{pvalue} - p-values
#'   estimated based on permutations and \code{adj_pvalue} - Benjamini &
#'   Hochberg adjusted p-values.
#'   
#' @examples 
#' 
#' #############################
#' ### sQTL analysis
#' #############################
#' # If possible, use BPPARAM = BiocParallel::MulticoreParam() with more workers
#' 
#' d <- data_dmSQTLdata
#' \donttest{
#' ### Filtering
#' d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_expr = 5, 
#'    min_samps_feature_prop = 0, minor_allele_freq = 5, 
#'    BPPARAM = BiocParallel::SerialParam())
#' 
#' ### Calculate dispersion
#' d <- dmDispersion(d, BPPARAM = BiocParallel::SerialParam())
#' 
#' ### Fit full model proportions
#' d <- dmFit(d, BPPARAM = BiocParallel::SerialParam())
#' 
#' ### Fit null model proportions and test for sQTLs
#' d <- dmTest(d, BPPARAM = BiocParallel::SerialParam())
#' plotPValues(d)
#' 
#' head(results(d))
#' 
#' }
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmSQTLdata}}, \code{\linkS4class{dmSQTLdata}}, 
#'   \code{\linkS4class{dmSQTLdispersion}}, \code{\linkS4class{dmSQTLfit}}
setClass("dmSQTLtest", 
  contains = "dmSQTLfit",
  representation(fit_null = "list",
    results = "data.frame"))


#####################################


setValidity("dmSQTLtest", function(object){
  # has to return TRUE when valid object!
  
  if(!length(object@counts) == length(object@fit_null))
    return("Different number of genes in 'counts' and 'fit_null'")
  
  if(!all(lapply(object@fit_null, class) == "MatrixList"))
    return("'fit_null' must be a list of MatrixLists")
  
  if(!nrow(object@results) == nrow(object@blocks))
    return("Different number of gene-SNP pairs in 'results' and in 'blocks'")
  
  return(TRUE)
  
})


################################################################################
### show and accessing methods
################################################################################


#' @rdname dmSQTLtest-class
#' @export
setMethod("results", "dmSQTLtest", function(x) x@results )



#################################

setMethod("show", "dmSQTLtest", function(object){
  
  callNextMethod(object)
  
  cat("* data accessors: results()\n")
  
})


###############################################################################
### dmTest
###############################################################################

#' @param permutations Character specifying which permutation scheme to apply for p-value calculation. When equal to \code{"all_genes"}, null distribution of p-values is calculated from all genes and the maximum number of permutation cycles is 10. When  \code{permutations = "per_gene"}, null distribution of p-values is calculated for each gene separately based on permutations of this individual gene. The latter approach may take a lot of computational time. We suggest using the first option.
#' @rdname dmTest
#' @export
setMethod("dmTest", "dmSQTLfit", function(x, permutations = "all_genes", 
  prop_mode = "constrOptim", 
  prop_tol = 1e-12, verbose = 0, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  stopifnot(length(prop_mode) == 1)
  stopifnot(prop_mode %in% c("constrOptim"))
  stopifnot(length(prop_tol) == 1)
  stopifnot(is.numeric(prop_tol) && prop_tol > 0)
  stopifnot(verbose %in% 0:2)
  stopifnot(permutations %in% c("all_genes", "per_gene"))
  
  
  fit_null <- dmSQTL_fitOneModel(counts = x@counts, genotypes = x@genotypes, 
    dispersion = slot(x, x@dispersion), model = "null", prop_mode = prop_mode, 
    prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  
  ### For "per_gene" approach, results is returned as a list, 
  ### so I do not have to split it in dmSQTL_permutations_per_gene
  if(permutations == "all_genes")
    return_list <- FALSE
  if(permutations == "per_gene")
    return_list <- TRUE
  
  
  results <- dmSQTL_test(fit_full = x@fit_full, fit_null = fit_null, 
    return_list = return_list, 
    verbose = verbose, BPPARAM = BPPARAM)
  
  
  if(verbose)
    message("\n** Running permutations..\n")
  
  ### Calculate adjusted p-values using permutations
  
  # P-value for a gene computed using all the permutations
  if(permutations == "all_genes")
    pval_adj_perm <- dmSQTL_permutations_all_genes(x = x, fit_null = fit_null, 
      results = results, max_nr_perm_cycles = 10, max_nr_min_nr_sign_pval = 1e3, 
      prop_mode = prop_mode, prop_tol = prop_tol, 
      verbose = verbose, BPPARAM = BPPARAM)
  
  # P-value for a gene computed using permutations of that gene 
  if(permutations == "per_gene")
    pval_adj_perm <- dmSQTL_permutations_per_gene(x = x, fit_null = fit_null, 
      results = results, max_nr_perm = 1e6, max_nr_sign_pval = 1e2, 
      prop_mode = prop_mode, prop_tol = prop_tol, 
      verbose = verbose, BPPARAM = BPPARAM)
  
  
  results$pvalue <- pval_adj_perm
  results$adj_pvalue <- p.adjust(pval_adj_perm, method="BH")
  
  
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





