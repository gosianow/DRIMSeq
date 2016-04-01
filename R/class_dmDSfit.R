#' @include class_dmDSdispersion.R
NULL

################################################################################
### dmDSfit class
################################################################################

#' dmDSfit object
#' 
#' dmDSfit extends the \code{\linkS4class{dmDSdispersion}} class by adding the
#' full model Dirichlet-multinomial feature proportion estimates needed for the
#' differential splicing analysis. Feature ratios are estimated for each gene
#' and each condition. Result of \code{\link{dmFit}}.
#' 
#' @return
#' 
#' \itemize{ \item \code{proportions(x)}: Get a data frame with estimated
#' feature ratios for each condition. \item \code{statistics(x)}: Get a data
#' frame with maximum log-likelihoods for each condition. }
#' 
#' @param x dmDSdispersion object.
#' @param ... Other parameters that can be defined by methods using this
#'   generic.
#'   
#' @slot dispersion Character specifying which type of dispersion was used for
#'   fitting: \code{"common_dispersion"} or \code{"genewise_dispersion"}.
#' @slot fit_full \code{\linkS4class{MatrixList}} containing the per gene
#'   feature ratios. Columns correspond to different conditions. Additionally,
#'   the full model likelihoods are stored in \code{metadata} slot.
#'   
#' @examples 
#' ###################################
#' ### Differential splicing analysis
#' ###################################
#' # If possible, use BPPARAM = BiocParallel::MulticoreParam() with more workers
#' 
#' d <- data_dmDSdata
#' \donttest{
#' ### Filtering
#' # Check what is the minimal number of replicates per condition 
#' table(samples(d)$group)
#' d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3, 
#'  min_samps_feature_prop = 0)
#' 
#' ### Calculate dispersion
#' d <- dmDispersion(d, BPPARAM = BiocParallel::SerialParam())
#' 
#' ### Fit full model proportions
#' d <- dmFit(d, BPPARAM = BiocParallel::SrialParam())
#' 
#' head(proportions(d))
#' head(statistics(d))
#' 
#' }
#' 
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmDSdata}}, \code{\linkS4class{dmDSdata}},
#'   \code{\linkS4class{dmDSdispersion}}, \code{\linkS4class{dmDStest}}
setClass("dmDSfit", 
  contains = "dmDSdispersion",
  representation(dispersion = "character",
    fit_full = "MatrixList"))


#####################################

setValidity("dmDSfit", function(object){
  # has to return TRUE when valid object!
  
  if(!length(object@dispersion) == 1)
    return("'dispersion' must have length 1")
  
  if(!object@dispersion %in% c("common_dispersion", "genewise_dispersion"))
    return("'dispersion' can have values 'common_dispersion' 
      or 'genewise_dispersion'")
  
  if(!length(object@counts) == length(object@fit_full))
    return("Different length of 'counts' and 'fit_full'")
  
  if(!ncol(object@fit_full) == nlevels(object@samples$group))
    return("Wrong number of groups in 'fit_full'")
  
  if(!ncol(object@fit_full@metadata) == nlevels(object@samples$group))
    return("Wrong number of groups in 'fit_full@metadata'")
  
  return(TRUE)
  
})



################################################################################
### accessing methods
################################################################################


#' @rdname dmDSfit-class
#' @export
setGeneric("proportions", function(x, ...) standardGeneric("proportions"))

#' @rdname dmDSfit-class
#' @export
setMethod("proportions", "dmDSfit", function(x){
  
  data.frame(gene_id = rep(names(x@counts), elementNROWS(x@counts)), 
    feature_id = rownames(x@counts@unlistData), x@fit_full@unlistData, 
    stringsAsFactors = FALSE, row.names = NULL)
  
})


#' @rdname dmDSfit-class
#' @export
setGeneric("statistics", function(x, ...) standardGeneric("statistics"))

#' @rdname dmDSfit-class
#' @export
setMethod("statistics", "dmDSfit", function(x){
  
  df <- data.frame(gene_id = names(x@counts), x@fit_full@metadata, 
    stringsAsFactors = FALSE, row.names = NULL)
  colnames(df)[-1] <- paste0("lik_", colnames(x@fit_full@metadata))
  return(df)
  
})



#################################

setMethod("show", "dmDSfit", function(object){
  
  callNextMethod(object)
  
  cat("  proportions(), statistics()\n")
  
  
})

################################################################################
### dmFit
################################################################################

#' Estimate proportions in Dirichlet-multinomial model
#' 
#' Maximum likelihood estimates of genomic feature (for instance, transcript,
#' exon, exonic bin) proportions in full Dirichlet-multinomial model used in
#' differential splicing or sQTL analysis. Full model estimation means that
#' proportions are estimated for every group/condition separately.
#' 
#' @param x \code{\linkS4class{dmDSdispersion}} or
#'   \code{\linkS4class{dmSQTLdispersion}} object.
#' @param ... Other parameters that can be defined by methods using this
#'   generic.
#' @export
setGeneric("dmFit", function(x, ...) standardGeneric("dmFit"))


####################################


#' @inheritParams dmDispersion
#' @param dispersion Character defining which dispersion should be used for 
#'   fitting. Possible values \code{"genewise_dispersion"} or 
#'   \code{"common_dispersion"}.
#' @return Returns a \code{\linkS4class{dmDSfit}} or 
#'   \code{\linkS4class{dmSQTLfit}} object.
#' @examples 
#' ###################################
#' ### Differential splicing analysis
#' ###################################
#' # If possible, use BPPARAM = BiocParallel::MulticoreParam() with more workers
#' 
#' d <- data_dmDSdata
#' \donttest{
#' ### Filtering
#' # Check what is the minimal number of replicates per condition 
#' table(samples(d)$group)
#' d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3, 
#'  min_samps_feature_prop = 0)
#' 
#' ### Calculate dispersion
#' d <- dmDispersion(d, BPPARAM = BiocParallel::SerialParam())
#' 
#' ### Fit full model proportions
#' d <- dmFit(d, BPPARAM = BiocParallel::SerialParam())
#' 
#' head(proportions(d))
#' head(statistics(d))
#' 
#' }
#' 
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmDSdata}}, \code{\link{data_dmSQTLdata}}, 
#'   \code{\link{plotFit}}, \code{\link{dmDispersion}}, \code{\link{dmTest}}
#' @rdname dmFit
#' @export
setMethod("dmFit", "dmDSdispersion", function(x, 
  dispersion = "genewise_dispersion", 
  prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = 0, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  stopifnot(length(dispersion) == 1)
  stopifnot(dispersion %in% c("genewise_dispersion", "common_dispersion"))
  stopifnot(length(prop_mode) == 1)
  stopifnot(prop_mode %in% c("constrOptimG", "constrOptim"))
  stopifnot(length(prop_tol) == 1)
  stopifnot(is.numeric(prop_tol) && prop_tol > 0)
  stopifnot(verbose %in% 0:2)
  
  fit_full <- dmDS_fitOneModel(counts = x@counts, samples = x@samples, 
    dispersion = slot(x, dispersion), model = "full", prop_mode = prop_mode, 
    prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
  
  return(new("dmDSfit", dispersion = dispersion, fit_full = fit_full,
    mean_expression = x@mean_expression, 
    common_dispersion = x@common_dispersion, 
    genewise_dispersion = x@genewise_dispersion, counts = x@counts, 
    samples = x@samples))
  
})


################################################################################
### plotFit
################################################################################

#' Plot feature proportions
#' 
#' @return Plot, per gene, the observed and estimated with Dirichlet-multinomial
#' model feature ratios. Estimated proportions are marked with diamond shapes.
#' 
#' @param x \code{\linkS4class{dmDSfit}}, \code{\linkS4class{dmDStest}} or
#'   \code{\linkS4class{dmSQTLfit}}, \code{\linkS4class{dmSQTLtest}} object.
#' @param ... Other parameters that can be defined by methods using this
#'   generic.
#' @export
setGeneric("plotFit", function(x, ...) standardGeneric("plotFit"))


#####################################


#' @inheritParams plotData
#' @param gene_id Character indicating a gene ID to be plotted.
#' @param plot_type Character defining the type of the plot produced. Possible
#'   values \code{"barplot"}, \code{"boxplot1"}, \code{"boxplot2"},
#'   \code{"lineplot"}, \code{"ribbonplot"}.
#' @param order Logical. Whether to plot the features ordered by their
#'   expression.
#' @param plot_full Logical. Whether to plot the proportions estimated by the
#'   full model.
#' @param plot_main Logical. Whether to plot a title with the information about
#'   the Dirichlet-multinomial estimates.
#'   
#' @examples 
#' 
#' ###################################
#' ### Differential splicing analysis
#' ###################################
#' # If possible, use BPPARAM = BiocParallel::MulticoreParam() with more workers
#' 
#' d <- data_dmDSdata
#' \donttest{
#' ### Filtering
#' # Check what is the minimal number of replicates per condition 
#' table(samples(d)$group)
#' d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3, 
#'  min_samps_feature_prop = 0)
#' 
#' ### Calculate dispersion
#' d <- dmDispersion(d, BPPARAM = BiocParallel::SerialParam())
#' 
#' ### Fit full model proportions
#' d <- dmFit(d, BPPARAM = BiocParallel::SerialParam())
#' 
#' ### Fit null model proportions and test for DS
#' d <- dmTest(d, BPPARAM = BiocParallel::SerialParam())
#' 
#' ### Plot feature proportions for top DS gene
#' res <- results(d)
#' res <- res[order(res$pvalue, decreasing = FALSE), ]
#' 
#' gene_id <- res$gene_id[1]
#' 
#' plotFit(d, gene_id = gene_id)
#' plotFit(d, gene_id = gene_id, plot_type = "lineplot")
#' plotFit(d, gene_id = gene_id, plot_type = "ribbonplot")
#' 
#' }
#' 
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmDSdata}}, \code{\link{data_dmSQTLdata}},
#'   \code{\link{plotData}}, \code{\link{plotDispersion}},
#'   \code{\link{plotTest}}
#' @rdname plotFit
#' @export
setMethod("plotFit", "dmDSfit", function(x, gene_id, plot_type = "barplot", 
  order = TRUE, plot_full = TRUE, plot_main = TRUE, out_dir = NULL){
  
  stopifnot(gene_id %in% names(x@counts))
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", 
    "ribbonplot"))
  stopifnot(is.logical(order))
  stopifnot(is.logical(plot_full))
  stopifnot(is.logical(plot_main))
  
  dmDS_plotFit(gene_id = gene_id, counts = x@counts, samples = x@samples, 
    dispersion = slot(x, x@dispersion), proportions_full = x@fit_full, 
    proportions_null = NULL, table = NULL, plot_type = plot_type, 
    order = order, plot_full = plot_full, plot_null = FALSE, 
    plot_main = plot_main, out_dir = out_dir)
  
  
})





