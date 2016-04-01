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
#'   statistics, \code{df} - degrees of freedom, \code{pvalue} - p-values and
#'   \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values.
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
#' plotTest(d)
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

#' @rdname dmTest
#' @export
setMethod("dmTest", "dmSQTLfit", function(x, prop_mode = "constrOptimG", 
  prop_tol = 1e-12, verbose = 0, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  stopifnot(length(prop_mode) == 1)
  stopifnot(prop_mode %in% c("constrOptimG", "constrOptim"))
  stopifnot(length(prop_tol) == 1)
  stopifnot(is.numeric(prop_tol) && prop_tol > 0)
  stopifnot(verbose %in% 0:2)
  
  fit_null <- dmSQTL_fitOneModel(counts = x@counts, genotypes = x@genotypes, 
    dispersion = slot(x, x@dispersion), model = "null", prop_mode = prop_mode, 
    prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  
  results <- dmSQTL_test(fit_full = x@fit_full, fit_null = fit_null, 
    verbose = verbose, BPPARAM = BPPARAM)
  
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
    res_new <- cbind(res[matching, c("gene_id", "block_id"), drop = FALSE], 
      snp_id, 
      res[matching, c("lr", "df", "pvalue", "adj_pvalue"), drop = FALSE])
    
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
### plotTest
###############################################################################

#' @rdname plotTest
#' @export
#' @importFrom grDevices pdf dev.off
setMethod("plotTest", "dmSQTLtest", function(x, out_dir = NULL){
  
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


################################################################################
### plotFit
################################################################################

#' @rdname plotFit
#' @export
setMethod("plotFit", "dmSQTLtest", function(x, gene_id, snp_id, 
  plot_type = "boxplot1", order = TRUE, plot_full = TRUE, plot_null = TRUE, 
  plot_main = TRUE, out_dir = NULL){
  
  stopifnot(gene_id %in% names(x@blocks))
  
  if(!snp_id %in% x@blocks[[gene_id, "snp_id"]])
    stop(paste0("gene ",gene_id, " and SNP ", snp_id, " do not match!"))
  
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", 
    "ribbonplot"))
  stopifnot(is.logical(order))
  stopifnot(is.logical(plot_full))
  stopifnot(is.logical(plot_null))
  stopifnot(is.logical(plot_main))
  
  
  dmSQTL_plotFit(gene_id = gene_id, snp_id = snp_id, counts = x@counts, 
    genotypes = x@genotypes, blocks = x@blocks, samples = x@samples, 
    dispersion = slot(x, x@dispersion), fit_full = x@fit_full, 
    fit_null = x@fit_null, table = x@results, plot_type = plot_type, 
    order = order, plot_full = plot_full, plot_null = plot_null, 
    plot_main = plot_main, out_dir = out_dir)
  
  
})





