#' @include class_dmSQTLdispersion.R class_dmDSfit.R
NULL

################################################################################
### dmSQTLfit class
################################################################################

#' dmSQTLfit object
#' 
#' dmSQTLfit extends the \code{\linkS4class{dmDSdispersion}} class by adding the
#' full model Dirichlet-multinomial feature proportion estimates needed for the
#' sQTL analysis. Feature ratios are estimated for each gene and each group that
#' is defined by different SNPs/blocks. Result of \code{\link{dmFit}}.
#' 
#' @slot dispersion Character specifying which type of dispersion was used for
#'   fitting: \code{"common_dispersion"} or \code{"genewise_dispersion"}.
#' @slot fit_full List of \code{\linkS4class{MatrixList}} objects. Each element
#'   of this list contains the full model proportion estimates for all the
#'   blocks associated with a given gene. Columns of MatrixLists correspond to 3
#'   genotypes (0,1,2). The full model likelihoods are stored in \code{metadata}
#'   slot.
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmSQTLdata}}, \code{\linkS4class{dmSQTLdata}},
#'   \code{\linkS4class{dmSQTLdispersion}}, \code{\linkS4class{dmSQTLtest}}
setClass("dmSQTLfit", 
  contains = "dmSQTLdispersion",
  representation(dispersion = "character",
    fit_full = "list"))

########################################


setValidity("dmSQTLfit", function(object){
  # has to return TRUE when valid object!
  
  if(!length(object@dispersion) == 1)
    return("'dispersion' must have length 1")
  
  if(!object@dispersion %in% c("common_dispersion", "genewise_dispersion"))
    return("'dispersion' can have values 'common_dispersion' or 'genewise_dispersion'")
  
  if(!length(object@counts) == length(object@fit_full))
    return("Different number of genes in 'counts' and 'fit_full'")
  
  if(!all(lapply(object@fit_full, class) == "MatrixList"))
    return("'fit_full' must be a list of MatrixLists")
  
  return(TRUE)
  
})


################################################################################
### show methods
################################################################################

setMethod("show", "dmSQTLfit", function(object){
  
  callNextMethod(object)
  
})


################################################################################
### dmFit
################################################################################


#' @rdname dmFit
#' @export
setMethod("dmFit", "dmSQTLdispersion", function(x, 
  dispersion = "genewise_dispersion", prop_mode = "constrOptim", 
  prop_tol = 1e-12, verbose = 0, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  stopifnot(length(dispersion) == 1)
  stopifnot(dispersion %in% c("genewise_dispersion", "common_dispersion"))
  stopifnot(length(prop_mode) == 1)
  stopifnot(prop_mode %in% c("constrOptim"))
  stopifnot(length(prop_tol) == 1)
  stopifnot(is.numeric(prop_tol) && prop_tol > 0)
  stopifnot(verbose %in% 0:2)
  
  fit_full <- dmSQTL_fitOneModel(counts = x@counts, genotypes = x@genotypes, 
    dispersion = slot(x, dispersion), model = "full", prop_mode = prop_mode, 
    prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  
  
  return(new("dmSQTLfit", dispersion = dispersion, fit_full = fit_full, 
    mean_expression = x@mean_expression, common_dispersion = x@common_dispersion, 
    genewise_dispersion = x@genewise_dispersion, counts = x@counts, 
    genotypes = x@genotypes, blocks = x@blocks, samples = x@samples))
  
  
})


################################################################################
### plotProportions
################################################################################


#' @param snp_id Character indicating a SNP ID to be plotted. \code{snp_id} must
#'   match \code{gene_id}.
#' @rdname plotProportions
#' @export
setMethod("plotProportions", "dmSQTLfit", function(x, gene_id, snp_id, 
  plot_type = "boxplot1", order = TRUE, plot_full = TRUE, 
  plot_main = TRUE, out_dir = NULL){
  
  stopifnot(gene_id %in% names(x@blocks))
  
  if(!snp_id %in% x@blocks[[gene_id, "snp_id"]])
    stop(paste0("gene ",gene_id, " and SNP ", snp_id, " do not match!"))
  
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", 
    "ribbonplot"))
  stopifnot(is.logical(order))
  stopifnot(is.logical(plot_full))
  stopifnot(is.logical(plot_main))
  
  dmSQTL_plotFit(gene_id = gene_id, snp_id = snp_id, counts = x@counts, 
    genotypes = x@genotypes, blocks = x@blocks, samples = x@samples, 
    dispersion = slot(x, x@dispersion), fit_full = x@fit_full, 
    fit_null = NULL, table = NULL, plot_type = plot_type, order = order, 
    plot_full = plot_full, plot_null = FALSE, plot_main = plot_main, 
    out_dir = out_dir)
  
  
})





