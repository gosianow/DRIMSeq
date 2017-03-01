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
#' feature ratios for each sample. \item \code{coefficients(x)}:  }
#' 
#' @param x dmDSdispersion object.
#' @param ... Other parameters that can be defined by methods using this 
#'   generic.
#'   
#' @slot design_fit_full Numeric matrix of the desing used to fit the full 
#'   model.
#' @slot fit_full \code{\linkS4class{MatrixList}} containing estimated feature 
#'   ratios in each sample based on the full Dirichlet-multinomial (DM) model.
#' @slot lik_full Numeric vector of the per gene DM full model likelihoods.
#' @slot coef_full \code{\linkS4class{MatrixList}} with the regression 
#'   coefficients based on the DM model
#' @slot fit_full_bb \code{\linkS4class{MatrixList}} containing estimated
#'   feature ratios in each sample based on the full beta-binomial (BB) model.
#' @slot lik_full_bb Numeric vector of the per gene BB full model likelihoods.
#' @slot coef_full_bb \code{\linkS4class{MatrixList}} with the regression 
#'   coefficients based on the BB model
#'   
#' @examples 
#' 
#' ###################################
#' ### Differential splicing analysis
#' ###################################
#' 
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmDSdata}}, \code{\linkS4class{dmDSdata}}, 
#'   \code{\linkS4class{dmDSdispersion}}, \code{\linkS4class{dmDStest}}
setClass("dmDSfit", 
  contains = "dmDSdispersion",
  representation(design_fit_full = "matrix",
    fit_full = "MatrixList",
    lik_full = "numeric",
    coef_full = "MatrixList",
    fit_full_bb = "MatrixList",
    lik_full_bb = "numeric",
    coef_full_bb = "MatrixList"))


# ------------------------------------------------------------------------------

setValidity("dmDSfit", function(object){
  # Has to return TRUE for a valid object!
  
  if(nrow(object@design_fit_full) == ncol(object@counts)){
    out <- TRUE
  }else{
    return(paste0("Number of rows in the design matrix must be equal 
      to the number of columns in counts"))
  }
  
  if(!length(object@fit_full) == length(object@counts))
    return("Different length of 'counts' and 'fit_full'")
  
  if(!length(object@lik_full) == length(object@counts))
    return("Different length of 'counts' and 'lik_full'")
  
  if(!length(object@coef_full) == length(object@counts))
    return("Different length of 'counts' and 'coef_full'")
  
  # TODO: Add more checks
  
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
  
  data.frame(gene_id = rep.int(names(x@fit_full), elementNROWS(x@fit_full)), 
    feature_id = rownames(x@fit_full@unlistData), x@fit_full@unlistData, 
    stringsAsFactors = FALSE, row.names = NULL)
  
})


# Generic for coefficients already exists in the stats package

#' @rdname dmDSfit-class
#' @param level Character specifying which type of results to return. Possible
#'   values \code{"gene"} or \code{"feature"}.
#' @export
setMethod("coefficients", "dmDSfit", function(object, level = "gene"){
  
  stopifnot(length(level) == 1)
  stopifnot(level %in% c("gene", "feature"))
  
  if(level == "gene"){
    out <- data.frame(gene_id = rep.int(names(object@coef_full), 
      elementNROWS(object@coef_full)), 
      feature_id = rownames(object@coef_full@unlistData), 
      object@coef_full@unlistData, 
      stringsAsFactors = FALSE, row.names = NULL)
  }
  if(level == "feature"){
    out <- data.frame(gene_id = rep.int(names(object@coef_full_bb), 
      elementNROWS(object@coef_full_bb)), 
      feature_id = rownames(object@coef_full_bb@unlistData), 
      object@coef_full_bb@unlistData, 
      stringsAsFactors = FALSE, row.names = NULL)
  }
  
  return(out)
  
})



# ------------------------------------------------------------------------------

setMethod("show", "dmDSfit", function(object){
  
  callNextMethod(object)
  
  cat("  proportions(), coefficients()\n")
  
  
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


# -----------------------------------------------------------------------------


#' @inheritParams dmDispersion
#' @param design Numeric matrix definig the full model.
#' @param bb_model Logical. Whether to perform the feature-level analysis using
#'   the beta-binomial model.
#' @return Returns a \code{\linkS4class{dmDSfit}} or 
#'   \code{\linkS4class{dmSQTLfit}} object.
#' @examples 
#' ###################################
#' ### Differential splicing analysis
#' ###################################
#' 
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmDSdata}}, \code{\link{data_dmSQTLdata}}, 
#'   \code{\link{plotProportions}}, \code{\link{dmDispersion}}, \code{\link{dmTest}}
#' @rdname dmFit
#' @export
setMethod("dmFit", "dmDSdispersion", function(x, design, 
  one_way = TRUE, bb_model = TRUE,
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = 0, BPPARAM = BiocParallel::SerialParam()){
  
  # Check design as in edgeR
  design <- as.matrix(design)
  stopifnot(nrow(design) == ncol(x@counts))
  
  ne <- limma::nonEstimable(design)
  if(!is.null(ne)) 
    stop(paste("Design matrix not of full rank. 
      The following coefficients not estimable:\n", paste(ne, collapse = " ")))
  
  if(!identical(x@design_dispersion, design))
    message(paste0("! The 'design' here is not identical as the 
      'design' used for dispersion estimation !\n"))
  
  # Check other parameters
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
  
  # Fit the DM model: proportions and likelihoods
  fit <- dmDS_fit(counts = x@counts, design = design, 
    dispersion = x@genewise_dispersion,
    one_way = one_way,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol,
    verbose = verbose, BPPARAM = BPPARAM)
  
  # Calculate the Beta-Binomial likelihoods for each feature
  if(bb_model){
    
    fit_bb <- bbDS_fit(counts = x@counts, fit = fit[["fit"]], design = design, 
      dispersion = x@genewise_dispersion,
      one_way = one_way,
      verbose = verbose, BPPARAM = BPPARAM)
    
    return(new("dmDSfit", design_fit_full = design, 
      fit_full = fit[["fit"]], lik_full = fit[["lik"]], coef_full = fit[["coef"]],
      lik_full_bb = fit_bb[["lik"]], coef_full_bb = fit_bb[["coef"]],
      mean_expression = x@mean_expression, 
      common_dispersion = x@common_dispersion, 
      genewise_dispersion = x@genewise_dispersion, 
      design_dispersion = x@design_dispersion,
      counts = x@counts, samples = x@samples))
    
  }else{
    
    return(new("dmDSfit", design_fit_full = design, 
      fit_full = fit[["fit"]], lik_full = fit[["lik"]], coef_full = fit[["coef"]],
      mean_expression = x@mean_expression, 
      common_dispersion = x@common_dispersion, 
      genewise_dispersion = x@genewise_dispersion, 
      design_dispersion = x@design_dispersion,
      counts = x@counts, samples = x@samples))
    
  }
  
  
})


################################################################################
### plotProportions
################################################################################

#' Plot feature proportions
#' 
#' This plot is available only for a group design, i.e., a comparison between 
#' groups.
#' 
#' @return Plot, per gene, the observed and estimated with Dirichlet-multinomial
#'   model feature proportions. Estimated group proportions are marked with
#'   diamond shapes.
#'   
#' @param x \code{\linkS4class{dmDSfit}}, \code{\linkS4class{dmDStest}} or 
#'   \code{\linkS4class{dmSQTLfit}}, \code{\linkS4class{dmSQTLtest}} object.
#' @param ... Other parameters that can be defined by methods using this 
#'   generic.
#' @export
setGeneric("plotProportions", function(x, ...) 
  standardGeneric("plotProportions"))


# ------------------------------------------------------------------------------


#' @inheritParams plotData
#' @param gene_id Character indicating a gene ID to be plotted.
#' @param group_variable Character indicating the grouping variable which is one
#'   of the columns in the \code{samples} slot of \code{x}.
#' @param plot_type Character defining the type of the plot produced. Possible 
#'   values \code{"barplot"}, \code{"boxplot1"}, \code{"boxplot2"}, 
#'   \code{"lineplot"}, \code{"ribbonplot"}.
#' @param order Logical. Whether to plot the features ordered by their 
#'   expression.
#' @param plot_full Logical. Whether to plot the proportions estimated by the 
#'   full model.
#' @param plot_main Logical. Whether to plot a title with the information about 
#'   the Dirichlet-multinomial estimates.
#' @param group_colors Character vector with colors for each group defined by 
#'   \code{group_variable}.
#' @param feature_colors Character vector with colors for each feature of gene
#'   defined by \code{gene_id}.
#'   
#' @examples 
#' 
#' ###################################
#' ### Differential splicing analysis
#' ###################################
#' 
#' 
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmDSdata}}, \code{\link{data_dmSQTLdata}}, 
#'   \code{\link{plotData}}, \code{\link{plotDispersion}}, 
#'   \code{\link{plotPValues}}
#' @rdname plotProportions
#' @export
setMethod("plotProportions", "dmDSfit", function(x, gene_id, group_variable, 
  plot_type = "barplot", order = TRUE, plot_fit = TRUE, plot_main = TRUE,
  group_colors = NULL, feature_colors = NULL){
  
  stopifnot(gene_id %in% names(x@counts))
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", 
    "ribbonplot"))
  stopifnot(is.logical(order))
  stopifnot(is.logical(plot_fit))
  stopifnot(is.logical(plot_main))
  stopifnot(length(group_variable) == 1)
  stopifnot(group_variable %in% colnames(samples(x)))

  group <- x@samples[, group_variable]
  counts_gene <- x@counts[[gene_id]]
  
  if(!is.null(group_colors) && 
      plot_type %in% c("barplot", "boxplot1", "lineplot"))
    stopifnot(length(group_colors) == nlevels(group))
  if(!is.null(feature_colors) && 
      plot_type %in% c("boxplot2", "ribbonplot"))
    stopifnot(length(feature_colors) == nrow(counts_gene))
  
  if(nrow(counts_gene) <= 1)
    stop("!Gene has to have at least 2 features! \n")
  
  # Order samples by group
  o <- order(group) 
  group <- group[o]
  counts_gene <- counts_gene[, o, drop = FALSE]
  
  main <- NULL
  
  if(plot_main){
    
    mean_expression_gene <- mean(colSums(counts_gene), na.rm = TRUE)
    
    main <- paste0(gene_id, "\n Mean expression = ", 
      round(mean_expression_gene))
    
    if(length(x@genewise_dispersion) > 0)
      dispersion_gene <- x@genewise_dispersion[gene_id]
    else
      dispersion_gene <- x@common_dispersion
    
    main <- paste0(main, ", Dispersion = ", round(dispersion_gene, 2))
    
  }
  
  
  prop_full <- NULL
  
  if(plot_fit){
    # Check if the design is equivalent to a oneway layout
    groups <- edgeR::designAsFactor(x@design_fit_full)
    
    if(nlevels(groups) == ncol(x@design_fit_full)){
      
      prop_full <- x@fit_full[[gene_id]][, !duplicated(group), drop = FALSE]
      colnames(prop_full) <- levels(group)
      
    }else{
      message("Fitted values are not plotted because the design does not 
        correspond to a group comparison defined by 'group_var'!")
    }
    
  }
  
  ggp <- dm_plotProportions(counts = counts_gene, group = group, 
    prop_full = prop_full, main = main, plot_type = plot_type, 
    order = order, group_colors = group_colors, feature_colors = feature_colors)
  
  return(ggp)  
  
})


























