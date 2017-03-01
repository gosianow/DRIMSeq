#' @include class_dmDSfit.R
NULL

###############################################################################
### dmDStest class
###############################################################################

#' dmDStest object
#' 
#' dmDStest extends the \code{\linkS4class{dmDSfit}} class by adding the null 
#' model Dirichlet-multinomial feature proportion estimates and the results of 
#' testing for differential splicing. Proportions are calculated for each gene 
#' from pooled (no grouping into conditions) counts. Result of 
#' \code{\link{dmTest}}.
#' 
#' @return
#' 
#' \itemize{ \item \code{proportions(x)}: Get a data frame with estimated 
#' feature ratios for full model and null models specified in 
#' \code{\link{dmTest}} with \code{compared_groups} parameter. \item 
#' \code{statistics(x)}: Get a data frame with full and null log-likelihoods and
#' degrees of freedom. \item \code{results(x)}: Get a data frame with results. 
#' See Slots. }
#' 
#' @param x dmDStest object.
#' @param ... Other parameters that can be defined by methods using this 
#'   generic.
#'   
#' @slot compared_groups Character vector specifying which groups/conditions 
#'   should be compared. By default, the comparison is done among all the groups
#'   specified by \code{group} column in \code{samples(x)}.
#' @slot fit_null \code{\linkS4class{MatrixList}}. Contains null proportions, 
#'   likelihoods and degrees of freedom for a comparison specified with 
#'   \code{compared_groups}.
#' @slot results Data frame with \code{gene_id} - gene IDs, \code{lr} - 
#'   likelihood ratio statistics, \code{df} - degrees of freedom, \code{pvalue} 
#'   - p-values and \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values 
#'   for comparison specified in \code{compared_groups}.
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
#' ### Fit full model proportions 
#' d <- dmFit(d, BPPARAM = BiocParallel::SerialParam())
#' 
#' ### Fit null model proportions and test for DS 
#' d <- dmTest(d, BPPARAM = BiocParallel::SerialParam()) 
#' plotTest(d)
#' 
#' head(proportions(d)) 
#' head(statistics(d)) 
#' head(results(d))
#' 
#' }
#' 
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmDSdata}}, \code{\linkS4class{dmDSdata}}, 
#'   \code{\linkS4class{dmDSdispersion}}, \code{\linkS4class{dmDSfit}}
setClass("dmDStest", 
  contains = "dmDSfit",
  representation(design_fit_null = "matrix",
    lik_null = "numeric",
    lik_null_bb = "numeric",
    results_gene = "data.frame",
    results_feature = "data.frame"))


##################################


setValidity("dmDStest", function(object){
  # Has to return TRUE when valid object!
  
  # TODO: set some validation!
  
  return(TRUE)
  
})

###############################################################################
### accessing methods
###############################################################################

#' @rdname dmDStest-class
#' @export
setMethod("proportions", "dmDStest", function(x){
  
  data.frame(gene_id = rep(names(x@counts), elementNROWS(x@counts)), 
    feature_id = rownames(x@counts@unlistData), x@fit_full@unlistData, 
    x@fit_null@unlistData, stringsAsFactors = FALSE, row.names = NULL)
  
})

###################################

#' @rdname dmDStest-class
#' @export
setMethod("statistics", "dmDStest", function(x){
  
  df <- data.frame(gene_id = names(x@counts), x@lik_full, x@lik_null, 
    stringsAsFactors = FALSE, row.names = NULL)
  
  return(df)
  
})

###################################

#' @rdname dmDStest-class
#' @export
setGeneric("results", function(x, ...) standardGeneric("results"))

#' @rdname dmDStest-class
#' @export
setMethod("results", "dmDStest", function(x, level = "gene") slot(x, paste0("results_", level)))



###################################

setMethod("show", "dmDStest", function(object){
  
  callNextMethod(object)
  
  cat("  results()\n")
  
})

###############################################################################
### dmTest
###############################################################################

#' Likelihood ratio test
#' 
#' First, estimate the null Dirichlet-multinomial model proportions, i.e.,
#' feature ratios are estimated based on pooled (no grouping into conditions)
#' counts. Use the likelihood ratio statistic to test for the difference between
#' feature proportions in different groups to identify the differentially
#' spliced genes (differential splicing analysis) or the sQTLs (sQTL analysis).
#' 
#' @param x \code{\linkS4class{dmDSfit}} or \code{\linkS4class{dmSQTLfit}}
#'   object.
#' @param ... Other parameters that can be defined by methods using this
#'   generic.
#' @export
setGeneric("dmTest", function(x, ...) standardGeneric("dmTest"))


##################################


#' @inheritParams dmFit
#' @param compared_groups Vector that defines which experimental conditions
#'   should be tested for differential splicing. By default, we test for a
#'   difference between any of the groups specified in \code{samples(x)$group}.
#'   Values in this vector should indicate levels or numbers of levels in
#'   \code{samples(x)$group}.
#'   
#' @return Returns a \code{\linkS4class{dmDStest}} or
#'   \code{\linkS4class{dmSQTLtest}} object.
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
#' plotTest(d)
#' 
#' head(proportions(d))
#' head(statistics(d))
#' head(results(d))
#' 
#' }
#' 
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmDSdata}}, \code{\link{data_dmSQTLdata}},
#'   \code{\link{plotTest}}, \code{\link{dmDispersion}}, \code{\link{dmFit}}
#' @rdname dmTest
#' @export
setMethod("dmTest", "dmDSfit", function(x, 
  coef = NULL, design = NULL, contrast = NULL,
  prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = 0, 
  BPPARAM = BiocParallel::SerialParam()){
  
  # Check parameters
  stopifnot(length(prop_mode) == 1)
  stopifnot(prop_mode %in% c("constrOptimG", "constrOptim"))
  stopifnot(length(prop_tol) == 1)
  stopifnot(is.numeric(prop_tol) && prop_tol > 0)
  stopifnot(verbose %in% 0:2)
  
  if(!sum(!unlist(lapply(list(coef, design, contrast), is.null))) == 1)
    stop(paste0("One of the ways to define the null model 'coef', 
      'design' or 'contrast' must be used!"))
  
  # Check coef
  if(!is.null(coef)){
    
    # Check the full model design matrix
    nbeta <- ncol(x@design_fit_full)
    if(nbeta < 2) 
      stop("Need at least two columns for design, usually the first is 
        the intercept column!")
    
    if(length(coef) > 1) 
      coef <- unique(coef)
    
    if(is.numeric(coef)){
      stopifnot(max(coef) <= nbeta)
    }else if(is.character(coef)){
      if(all(coef %in% colnames(x@design_fit_full))) 
        stop("'coef' does not match the columns of the design matrix!")
      coef <- match(coef, colnames(x@design_fit_full))
    }
    
    # Null design matrix
    design0 <- x@design_fit_full[, -coef, drop = FALSE]
    
  }
  
  # Check design
  if(!is.null(design)){
    
    # Check design as in edgeR
    design <- as.matrix(design)
    stopifnot(nrow(design) == ncol(x@counts))
    
    ne <- limma::nonEstimable(design)
    if(!is.null(ne)) 
      stop(paste("Design matrix not of full rank. 
        The following coefficients not estimable:\n", 
        paste(ne, collapse = " ")))
    
    # Null design matrix
    design0 <- design
    
  }
  
  # Check contrast exactly as in edgeR in glmLRT()
  if(!is.null(contrast)){
    
    contrast <- as.matrix(contrast)
    qrc <- qr(contrast)
    ncontrasts <- qrc$rank
    
    if(ncontrasts == 0) 
      stop("Contrasts are all zero!")
    
    coef <- 1:ncontrasts
    
    Dvec <- rep.int(1, nlibs)
    Dvec[coef] <- diag(qrc$qr)[coef]
    Q <- qr.Q(qrc, complete = TRUE, Dvec = Dvec)
    design <- design %*% Q
    
    # Null design matrix
    design0 <- design[, -coef, drop = FALSE]
    
  }
  
  # Fit the DM null model: proportions and likelihoods
  fit0 <- dmDS_fit(counts = x@counts, design = design0, 
    dispersion = x@genewise_dispersion,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    verbose = verbose, BPPARAM = BPPARAM)
  
  # Calculate the DM degrees of freedom for the LR test: df_full - df_null
  df <- (ncol(x@design_fit_full) - ncol(design0)) * 
    (elementNROWS(x@coef_full) - 1)
  
  results_gene <- dm_LRT(lik_full = x@lik_full, 
    lik_null = fit0[["lik"]], df = df, verbose = verbose)
  
  results_gene <- data.frame(gene_id = rownames(results_gene), 
    results_gene, stringsAsFactors = FALSE, row.names = NULL)
  
  
  # Calculate the Beta-Binomial null likelihoods for each feature
  fit0_bb <- bbDS_fit(counts = x@counts, fit = fit0[["fit"]], design = design0, 
    dispersion = x@genewise_dispersion,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    verbose = verbose, BPPARAM = BPPARAM)
  
  # Calculate the BB degrees of freedom for the LR test
  df <- rep(ncol(x@design_fit_full) - ncol(design0), length(x@lik_full_bb))
  
  results_feature <- dm_LRT(lik_full = x@lik_full_bb, 
    lik_null = fit0_bb[["lik"]], df = df, verbose = verbose)
  
  results_feature <- data.frame(
    gene_id = rep(names(x@counts), elementNROWS(x@counts)), 
    feature_id = rownames(results_feature), 
    results_feature, stringsAsFactors = FALSE, row.names = NULL)
  
  return(new("dmDStest", 
    results_gene = results_gene, results_feature = results_feature,
    design_fit_null = design0, 
    lik_null = fit0[["lik"]],
    lik_null_bb = fit0_bb[["lik"]],
    design_fit_full = x@design_fit_full, 
    fit_full = x@fit_full, lik_full = x@lik_full, coef_full = x@coef_full,
    lik_full_bb = x@lik_full_bb,  coef_full_bb = x@coef_full_bb,
    mean_expression = x@mean_expression, 
    common_dispersion = x@common_dispersion, 
    genewise_dispersion = x@genewise_dispersion, 
    design_dispersion = x@design_dispersion,
    counts = x@counts, samples = x@samples))
  
})


###############################################################################
### plotTest
###############################################################################

#' Plot p-values distribution
#' 
#' @return Plot a histogram of p-values.
#' 
#' @param x \code{\linkS4class{dmDStest}} or \code{\linkS4class{dmSQTLtest}}
#'   object.
#' @export
setGeneric("plotTest", function(x, ...) standardGeneric("plotTest"))



####################################

#' @inheritParams plotData
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
#' plotTest(d)
#' 
#' 
#' }
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmDSdata}}, \code{\link{data_dmSQTLdata}},
#'   \code{\link{plotData}}, \code{\link{plotDispersion}}, \code{\link{plotFit}}
#' @rdname plotTest
#' @export
#' @importFrom grDevices pdf dev.off
setMethod("plotTest", "dmDStest", function(x, out_dir = NULL){
  
  ggp <- dm_plotPvalues(pvalues = x@results[, "pvalue"])
  
  if(!is.null(out_dir)){
    pdf(paste0(out_dir, "hist_pvalues.pdf"))
    print(ggp)
    dev.off()
  }else{
    return(ggp)  
  }
  
})


###############################################################################
### plotFit
###############################################################################

#' @param plot_null Logical. Whether to plot the proportions estimated by the
#'   null model.
#' @rdname plotFit
#' @export
setMethod("plotFit", "dmDStest", function(x, gene_id, plot_type = "barplot", 
  order = TRUE, plot_full = TRUE, plot_null = TRUE, plot_main = TRUE, 
  out_dir = NULL){
  
  stopifnot(all(gene_id %in% names(x@counts)))
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", 
    "ribbonplot"))
  stopifnot(is.logical(order))
  stopifnot(is.logical(plot_full))
  stopifnot(is.logical(plot_null))
  stopifnot(is.logical(plot_main))
  
  
  compared_groups <- x@compared_groups
  
  samps <- x@samples$group %in% compared_groups
  
  samples = x@samples[samps, , drop = FALSE]
  samples$sample_id <- factor(samples$sample_id)
  samples$group <- factor(samples$group)
  
  results <- x@results
  
  dmDS_plotFit(gene_id = gene_id, counts = x@counts[, samps, drop = FALSE], 
    samples = samples, 
    dispersion = slot(x, x@dispersion), 
    proportions_full = x@fit_full[, compared_groups, drop = FALSE], 
    proportions_null = x@fit_null, table = results, plot_type = plot_type, 
    order = order, plot_full = plot_full, plot_null = plot_null, 
    plot_main = plot_main, out_dir = out_dir)
  
  
})



###############################################################################
### dmTwoStageTest
###############################################################################

#' Two-stage test
#' 
#' Two-stage test
#' 
#' @param x \code{\linkS4class{dmDStest}} or \code{\linkS4class{dmSQTLtest}}
#'   object.
#' @param ... Other parameters that can be defined by methods using this
#'   generic.
#' @export
setGeneric("dmTwoStageTest", function(x, ...) standardGeneric("dmTwoStageTest"))


#' @inheritParams dmTest
#' @param FDR Numeric. Cutoff for the FDR.
#' @return Returns a data frame with adjusted feature level p-values.
#' @author Malgorzata Nowicka
#' @rdname dmTwoStageTest
#' @export
setMethod("dmTwoStageTest", "dmDStest", function(x, FDR = 0.05, verbose = 0){
  
  stopifnot(verbose %in% 0:2)
  stopifnot(length(FDR) == 1)
  stopifnot(class(FDR) == "numeric")
  
  table <- dm_twoStageTest(pvalue_gene = x@results_gene, 
    pvalue_feature = x@results_feature, FDR = FDR, verbose = verbose)
  
  return(table)
  
})







































