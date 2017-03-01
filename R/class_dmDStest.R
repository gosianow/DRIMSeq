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
#' \itemize{ \item \code{results(x)}: Get a data frame with gene-level or
#' feature-level results. See Slots. }
#' 
#' @param x dmDStest object.
#' @param ... Other parameters that can be defined by methods using this 
#'   generic.
#'   
#' @slot design_fit_null Numeric matrix of the desing used to fit the null 
#'   model.
#' @slot lik_null Numeric vector of the per gene DM null model likelihoods.
#' @slot lik_null_bb Numeric vector of the per gene BB null model likelihoods.
#' @slot results_gene Data frame with the gene-level results including: 
#'   \code{gene_id} - gene IDs, \code{lr} - likelihood ratio statistics based on
#'   the DM model, \code{df} - degrees of freedom, \code{pvalue} - p-values and 
#'   \code{adj_pvalue} - Benjamini & Hochberg adjusted p-values.
#' @slot results_feature Data frame with the feature-level results including: 
#'   \code{gene_id} - gene IDs, \code{feature_id} - feature IDs, \code{lr} - 
#'   likelihood ratio statistics based on the BB model, \code{df} - degrees of 
#'   freedom, \code{pvalue} - p-values and \code{adj_pvalue} - Benjamini & 
#'   Hochberg adjusted p-values.
#'   
#' @examples 
#' 
#' ###################################
#' ### Differential splicing analysis
#' ###################################
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
  
  # TODO: Add more checks
  
  return(TRUE)
  
})

###############################################################################
### accessing methods
###############################################################################


#' @rdname dmDStest-class
#' @export
setGeneric("results", function(x, ...) standardGeneric("results"))

#' @rdname dmDStest-class
#' @param level Character specifying which type of results to return. Possible
#'   values \code{"gene"} or \code{"feature"}.
#' @export
setMethod("results", "dmDStest", function(x, level = "gene"){
  stopifnot(length(level) == 1)
  stopifnot(level %in% c("gene", "feature"))
  slot(x, paste0("results_", level))
})


# -----------------------------------------------------------------------------

setMethod("show", "dmDStest", function(object){
  
  callNextMethod(object)
  
  cat("  results()\n")
  
})

###############################################################################
### dmTest
###############################################################################

#' Likelihood ratio test
#' 
#' First, estimate the null Dirichlet-multinomial and beta-binomial model
#' parameters and likelihoods. Second, perform the gene-level (DM model) and
#' feature-level (BB model) likelihood ratio tests.
#' 
#' @param x \code{\linkS4class{dmDSfit}} or \code{\linkS4class{dmSQTLfit}} 
#'   object.
#' @param ... Other parameters that can be defined by methods using this 
#'   generic.
#' @export
setGeneric("dmTest", function(x, ...) standardGeneric("dmTest"))


##################################


#' @inheritParams dmFit
#' @param coef Integer or character vector indicating which coefficients of the 
#'   linear model are to be tested equal to zero. Values must indicate column 
#'   numbers or column names of the \code{design} used in 
#'   \code{\linkS4class{dmFit}}.
#' @param design Numeric matrix definig the null model.
#' @param contrast Numeric vector or matrix specifying one or more contrasts of 
#'   the linear model coefficients to be tested equal to zero. Number of rows 
#'   must equal to the number of columns of \code{design} used in 
#'   \code{\linkS4class{dmFit}}.
#'   
#' @details One must specify one of the arguments: \code{coef}, \code{design} or
#'   \code{contrast}.
#'   
#' @return Returns a \code{\linkS4class{dmDStest}} or 
#'   \code{\linkS4class{dmSQTLtest}} object.
#' @examples 
#' 
#' ###################################
#' ### Differential splicing analysis
#' ###################################
#' 
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmDSdata}}, \code{\link{data_dmSQTLdata}}, 
#'   \code{\link{plotPValues}}, \code{\link{dmDispersion}}, \code{\link{dmFit}}
#' @rdname dmTest
#' @export
setMethod("dmTest", "dmDSfit", function(x, 
  coef = NULL, design = NULL, contrast = NULL, 
  one_way = TRUE,
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
  stopifnot(coef_mode %in% c("optim", "nlminb", "Rcgmin"))
  stopifnot(length(coef_tol) == 1)
  stopifnot(is.numeric(coef_tol) && coef_tol > 0)
  
  stopifnot(verbose %in% 0:2)
  
  if(!sum(!unlist(lapply(list(coef, design, contrast), is.null))) == 1)
    stop(paste0("Only one of the ways to define the null model 'coef', 
      'design' or 'contrast' can be used!"))
  
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
    
    design <- x@design_fit_full 
    contrast <- as.matrix(contrast)
    stopifnot(nrow(contrast) == ncol(design))
    
    qrc <- qr(contrast)
    ncontrasts <- qrc$rank
    
    if(ncontrasts == 0) 
      stop("Contrasts are all zero!")
    
    coef <- 1:ncontrasts
    
    nlibs <- nrow(design)
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
    one_way = one_way,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol, 
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
    one_way = one_way,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol,
    verbose = verbose, BPPARAM = BPPARAM)
  
  # Calculate the BB degrees of freedom for the LR test
  df <- rep.int(ncol(x@design_fit_full) - ncol(design0), length(x@lik_full_bb))
  
  results_feature <- dm_LRT(lik_full = x@lik_full_bb, 
    lik_null = fit0_bb[["lik"]], df = df, verbose = verbose)
  
  results_feature <- data.frame(
    gene_id = rep.int(names(x@counts), elementNROWS(x@counts)), 
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
### plotPValues
###############################################################################

#' Plot p-values distribution
#' 
#' @return Plot a histogram of p-values.
#' 
#' @param x \code{\linkS4class{dmDStest}} or \code{\linkS4class{dmSQTLtest}}
#'   object.
#' @export
setGeneric("plotPValues", function(x, ...) standardGeneric("plotPValues"))



####################################

#' @inheritParams plotData
#' @examples
#' 
#' ###################################
#' ### Differential splicing analysis
#' ###################################
#' 
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmDSdata}}, \code{\link{data_dmSQTLdata}},
#'   \code{\link{plotData}}, \code{\link{plotDispersion}}, \code{\link{plotFit}}
#' @rdname plotPValues
#' @export
#' @importFrom grDevices pdf dev.off
setMethod("plotPValues", "dmDStest", function(x, out_dir = NULL){
  
  ggp <- dm_plotPValues(pvalues = x@results[, "pvalue"])
  
  if(!is.null(out_dir)){
    pdf(paste0(out_dir, "hist_pvalues.pdf"))
    print(ggp)
    dev.off()
  }else{
    return(ggp)  
  }
  
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
#' @return Returns a data frame with adjusted feature-level p-values.
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







































