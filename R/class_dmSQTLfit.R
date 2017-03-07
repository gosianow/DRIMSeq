#' @include class_dmSQTLprecision.R class_dmDSfit.R
NULL

################################################################################
### dmSQTLfit class
################################################################################

#' dmSQTLfit object
#' 
#' dmSQTLfit extends the \code{\linkS4class{dmSQTLprecision}} class by adding 
#' the full model Dirichlet-multinomial (DM) likelihoods,
#' regression coefficients and feature proportion estimates needed for the 
#' transcript/exon usage QTL analysis. Full model is defined by the genotype of
#' a SNP associated with a gene. Estimation takes place for all the genes and 
#' all the SNPs/blocks assigned to the genes. Result of \code{\link{dmFit}}.
#' 
#' @slot fit_full List of \code{\linkS4class{MatrixList}} objects containing 
#'   estimated feature ratios in each sample based on the full 
#'   Dirichlet-multinomial (DM) model.
#' @slot lik_full List of numeric vectors of the per gene DM full model 
#'   likelihoods.
#' @slot coef_full \code{\linkS4class{MatrixList}} with the regression 
#'   coefficients based on the DM model.
#' @examples 
#' # --------------------------------------------------------------------------
#' # Create dmSQTLdata object
#' # --------------------------------------------------------------------------
#' # Use subsets of data defined in the GeuvadisTranscriptExpr package
#' 
#' library(GeuvadisTranscriptExpr)
#' \donttest{
#' geuv_counts <- GeuvadisTranscriptExpr::counts
#' geuv_genotypes <- GeuvadisTranscriptExpr::genotypes
#' geuv_gene_ranges <- GeuvadisTranscriptExpr::gene_ranges
#' geuv_snp_ranges <- GeuvadisTranscriptExpr::snp_ranges
#' 
#' colnames(geuv_counts)[c(1,2)] <- c("feature_id", "gene_id")
#' colnames(geuv_genotypes)[4] <- "snp_id"
#' geuv_samples <- data.frame(sample_id = colnames(geuv_counts)[-c(1,2)])
#' 
#' d <- dmSQTLdata(counts = geuv_counts, gene_ranges = geuv_gene_ranges,  
#'   genotypes = geuv_genotypes, snp_ranges = geuv_snp_ranges, 
#'   samples = geuv_samples, window = 5e3)
#' 
#' # --------------------------------------------------------------------------
#' # sQTL analysis - simple group comparison
#' # --------------------------------------------------------------------------
#' 
#' ## Filtering
#' d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_expr = 5,
#'   minor_allele_freq = 5, min_gene_expr = 10, min_feature_expr = 10)
#'   
#' plotData(d)
#' 
#' ## To make the analysis reproducible
#' set.seed(123)
#' ## Calculate precision
#' d <- dmPrecision(d)
#' 
#' plotPrecision(d)
#' 
#' ## Fit full model proportions
#' d <- dmFit(d)
#' }
#' @author Malgorzata Nowicka
#' @seealso \code{\linkS4class{dmSQTLdata}}, 
#'   \code{\linkS4class{dmSQTLprecision}}, \code{\linkS4class{dmSQTLtest}}
setClass("dmSQTLfit", 
  contains = "dmSQTLprecision",
  representation(fit_full = "list",
    lik_full = "list",
    coef_full = "list"))

########################################


setValidity("dmSQTLfit", function(object){
  # Has to return TRUE when valid object
  
  # TODO: Add checks for other slots
  
  if(!length(object@counts) == length(object@lik_full))
    return("Different number of genes in 'counts' and 'lik_full'")
  
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


#' @details In the QTL analysis, currently, genotypes are defined as numeric
#' values 0, 1, and 2. When \code{one_way = TRUE}, simple multiple group fitting
#' is performed. When \code{one_way = FALSE}, a regression framework is used
#' with the design matrix defined by a formula \code{~ group} where group is a 
#' continuous (not categorical) variable with values 0, 1, and 2.
#' @rdname dmFit
#' @export
setMethod("dmFit", "dmSQTLprecision", function(x, one_way = TRUE, 
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = 0, BPPARAM = BiocParallel::SerialParam()){
  
  # Check parameters
  stopifnot(is.logical(one_way))
  
  stopifnot(length(prop_mode) == 1)
  stopifnot(prop_mode %in% c("constrOptim"))
  stopifnot(length(prop_tol) == 1)
  stopifnot(is.numeric(prop_tol) && prop_tol > 0)
  
  stopifnot(length(coef_mode) == 1)
  stopifnot(coef_mode %in% c("optim", "nlminb", "nlm"))
  stopifnot(length(coef_tol) == 1)
  stopifnot(is.numeric(coef_tol) && coef_tol > 0)
  
  stopifnot(verbose %in% 0:3)
  
  fit <- dmSQTL_fit(counts = x@counts, genotypes = x@genotypes, 
    precision = x@genewise_precision,
    one_way = one_way, group_formula = ~ group,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol,
    return_fit = FALSE, return_coef = FALSE,
    verbose = verbose, BPPARAM = BPPARAM)
  
  return(new("dmSQTLfit", lik_full = fit[["lik"]], fit_full = fit[["fit"]],
    mean_expression = x@mean_expression, 
    common_precision = x@common_precision, 
    genewise_precision = x@genewise_precision, 
    counts = x@counts, genotypes = x@genotypes,
    blocks = x@blocks, samples = x@samples))
  
})


################################################################################
### plotProportions
################################################################################


#' @param snp_id Character indicating the ID of a SNP to be plotted.
#' @details In the QTL analysis, plotting of fitted proportions is deactivated 
#'   even when \code{plot_fit = TRUE}. It is due to the fact that neither fitted
#'   values nor regression coefficients are returned by the \code{dmFit}
#'   function as they occupy a lot of memory.
#' @rdname plotProportions
#' @export
setMethod("plotProportions", "dmSQTLfit", function(x, gene_id, snp_id, 
  plot_type = "boxplot1", order_features = TRUE, order_samples = TRUE,
  plot_fit = FALSE, plot_main = TRUE, 
  group_colors = NULL, feature_colors = NULL){
  
  stopifnot(gene_id %in% names(x@blocks))
  
  if(!snp_id %in% x@blocks[[gene_id, "snp_id"]])
    stop(paste0("gene ",gene_id, " and SNP ", snp_id, " do not match!"))
  
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", 
    "ribbonplot"))
  stopifnot(is.logical(order_features))
  stopifnot(is.logical(order_samples))
  stopifnot(is.logical(plot_fit))
  stopifnot(is.logical(plot_main))
  
  counts_gene <- x@counts[[gene_id]]
  block_id <- x@blocks[[gene_id]][x@blocks[[gene_id]][, "snp_id"] == snp_id,
    "block_id"]
  group <- x@genotypes[[gene_id]][block_id, ] 
  
  if(!is.null(group_colors) && 
      plot_type %in% c("barplot", "boxplot1", "lineplot"))
    stopifnot(length(group_colors) == nlevels(group))
  if(!is.null(feature_colors) && 
      plot_type %in% c("boxplot2", "ribbonplot"))
    stopifnot(length(feature_colors) == nrow(counts_gene))
  
  if(nrow(counts_gene) <= 1)
    stop("!Gene has to have at least 2 features! \n")
  
  # Remove NAs
  nonNAs <- !(is.na(counts_gene[1,]) | is.na(group))
  counts_gene <- counts_gene[, nonNAs, drop = FALSE]
  group <- factor(group[nonNAs])
  
  main <- NULL
  
  if(plot_main){
    
    mean_expression_gene <- mean(colSums(counts_gene), na.rm = TRUE)
    
    main <- paste0(gene_id, " : ", snp_id, " : ", block_id,
      "\n Mean expression = ", round(mean_expression_gene))
    
    precision_gene <- x@genewise_precision[[gene_id]][block_id]
    
    main <- paste0(main, ", Precision = ", round(precision_gene, 2))
    
  }
  
  fit_full <- NULL
  
  if(plot_fit && length(x@fit_full) > 0){

    fit_full <- x@fit_full[[gene_id]][[which(rownames(x@genotypes[[gene_id]]) == 
        block_id)]][, nonNAs, drop = FALSE]

  }

  ggp <- dm_plotProportions(counts = counts_gene, group = group, 
    fit_full = fit_full, main = main, plot_type = plot_type, 
    order_features = order_features, order_samples = order_samples,
    group_colors = group_colors, feature_colors = feature_colors)
  
  return(ggp)  
  
})





