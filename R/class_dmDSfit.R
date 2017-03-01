#' @include class_dmDSprecision.R
NULL

################################################################################
### dmDSfit class
################################################################################

#' dmDSfit object
#' 
#' dmDSfit extends the \code{\linkS4class{dmDSprecision}} class by adding the 
#' full model Dirichlet-multinomial (DM) and beta-binomial (BB) likelihoods, 
#' regression coefficients and feature proportion estimates. Result of calling 
#' the \code{\link{dmFit}} function.
#' 
#' @return
#' 
#' \itemize{ \item \code{design(object)}: Get a matrix with the full design.
#' \item \code{proportions(x)}: Get a data frame with estimated feature ratios
#' for each sample. \item \code{coefficients(x)}: Get the DM or BB regression
#' coefficients. }
#' 
#' @param x,object dmDSprecision object.
#' @param ... Other parameters that can be defined by methods using this 
#'   generic.
#'   
#' @slot design_fit_full Numeric matrix of the design used to fit the full 
#'   model.
#' @slot fit_full \code{\linkS4class{MatrixList}} containing estimated feature 
#'   ratios in each sample based on the full Dirichlet-multinomial (DM) model.
#' @slot lik_full Numeric vector of the per gene DM full model likelihoods.
#' @slot coef_full \code{\linkS4class{MatrixList}} with the regression 
#'   coefficients based on the DM model.
#' @slot fit_full_bb \code{\linkS4class{MatrixList}} containing estimated 
#'   feature ratios in each sample based on the full beta-binomial (BB) model.
#' @slot lik_full_bb Numeric vector of the per gene BB full model likelihoods.
#' @slot coef_full_bb \code{\linkS4class{MatrixList}} with the regression 
#'   coefficients based on the BB model.
#'   
#' @examples 
#' # --------------------------------------------------------------------------
#' # Create dmDSdata object 
#' # --------------------------------------------------------------------------
#' ## Get kallisto transcript counts from the 'PasillaTranscriptExpr' package
#' 
#' library(PasillaTranscriptExpr)
#' \donttest{
#' data_dir  <- system.file("extdata", package = "PasillaTranscriptExpr")
#' 
#' ## Load metadata
#' pasilla_metadata <- read.table(file.path(data_dir, "metadata.txt"), 
#' header = TRUE, as.is = TRUE)
#' 
#' ## Load counts
#' pasilla_counts <- read.table(file.path(data_dir, "counts.txt"), 
#' header = TRUE, as.is = TRUE)
#' 
#' ## Create a pasilla_samples data frame
#' pasilla_samples <- data.frame(sample_id = pasilla_metadata$SampleName, 
#'   group = pasilla_metadata$condition)
#' levels(pasilla_samples$group)
#' 
#' ## Create a dmDSdata object
#' d <- dmDSdata(counts = pasilla_counts, samples = pasilla_samples)
#' 
#' ## Use a subset of genes, which is defined in the following file
#' gene_id_subset <- readLines(file.path(data_dir, "gene_id_subset.txt"))
#' 
#' d <- d[names(d) %in% gene_id_subset, ]
#' 
#' # --------------------------------------------------------------------------
#' # Differential transcript usage analysis - simple two group comparison 
#' # --------------------------------------------------------------------------
#' 
#' ## Filtering
#' ## Check what is the minimal number of replicates per condition 
#' table(samples(d)$group)
#' 
#' d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3,
#'   min_gene_expr = 10, min_feature_expr = 10)
#' 
#' plotData(d)
#' 
#' ## Create the design matrix
#' design_full <- model.matrix(~ group, data = samples(d))
#' 
#' ## To make the analysis reproducible
#' set.seed(123)
#' ## Calculate precision
#' d <- dmPrecision(d, design = design_full)
#' 
#' plotPrecision(d)
#' 
#' head(mean_expression(d))
#' common_precision(d)
#' head(genewise_precision(d))
#' 
#' ## Fit full model proportions
#' d <- dmFit(d, design = design_full)
#' 
#' ## Get fitted proportions
#' head(proportions(d))
#' ## Get the DM regression coefficients (gene-level) 
#' head(coefficients(d))
#' ## Get the BB regression coefficients (feature-level) 
#' head(coefficients(d), level = "feature")
#' }
#' @author Malgorzata Nowicka
#' @seealso \code{\linkS4class{dmDSdata}}, \code{\linkS4class{dmDSprecision}}, 
#'   \code{\linkS4class{dmDStest}}
setClass("dmDSfit", 
  contains = "dmDSprecision",
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
    return("Different number of genes in 'counts' and 'fit_full'")
  
  if(!length(object@lik_full) == length(object@counts))
    return("Different number of genes in 'counts' and 'lik_full'")
  
  if(!length(object@coef_full) == length(object@counts))
    return("Different number of genes in 'counts' and 'coef_full'")
  
  # TODO: Add more checks for BB
  
  return(TRUE)
  
})



################################################################################
### accessing methods
################################################################################

#' @rdname dmDSfit-class
#' @inheritParams dmDSprecision-class
#' @export
setMethod("design", "dmDSfit", function(object, type = "full_model"){
  
  stopifnot(type %in% c("precision", "full_model", "null_model"))
  
  if(type == "precision")
    object@design_precision
  else if(type == "full_model")
    object@design_fit_full
  else
    NULL
  
})


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

#' Fit the Dirichlet-multinomial and/or the beta-binomial full model regression
#' 
#' Obtain the maximum likelihood estimates of Dirichlet-multinomial (gene-level)
#' and/or beta-binomial (feature-level) regression coefficients, feature 
#' proportions in each sample and corresponding likelihoods. In the differential
#' exon/transcript usage analysis, the regression model is defined by a design 
#' matrix. In the exon/transcript usage QTL analysis, regression models are 
#' defined by genotypes. Currently, beta-binomial model is implemented only in
#' the differential usage analysis.
#' 
#' @param x \code{\linkS4class{dmDSprecision}} or 
#'   \code{\linkS4class{dmSQTLprecision}} object.
#' @param ... Other parameters that can be defined by methods using this 
#'   generic.
#' @export
setGeneric("dmFit", function(x, ...) standardGeneric("dmFit"))


# -----------------------------------------------------------------------------


#' @inheritParams dmPrecision
#'   
#' @details In the regression framework here, we adapt the idea from 
#'   \code{\link[edgeR]{glmFit}} in \code{\link{edgeR}} about using a shortcut 
#'   algorithm when the design is equivalent to simple group fitting. In such a 
#'   case, we estimate the DM proportions for each group of samples separately 
#'   and then recalculate the DM (and/or the BB) regression coefficients 
#'   corresponding to the design matrix. If the design matrix does not define a 
#'   simple group fitting, for example, when it contains a column with
#'   continuous values, then the regression framework is used to directly
#'   estimate the regression coefficients.
#'   
#'   Arguments that are used for the proportion estimation in each group when 
#'   the shortcut fitting can be used start with \code{prop_}, and those that 
#'   are used in the regression framework start with \code{coef_}.
#'   
#'   In the differential transcript usage analysis, setting \code{one_way = 
#'   TRUE} allows switching to the shortcut algorithm only if the design is 
#'   equivalent to simple group fitting. \code{one_way = FALSE} forces usage of 
#'   the regression framework.
#'   
#'   
#' @param design Numeric matrix defining the full model.
#' @param bb_model Logical. Whether to perform the feature-level analysis using 
#'   the beta-binomial model.
#' @return Returns a \code{\linkS4class{dmDSfit}} or 
#'   \code{\linkS4class{dmSQTLfit}} object.
#' @examples 
#' # --------------------------------------------------------------------------
#' # Create dmDSdata object 
#' # --------------------------------------------------------------------------
#' ## Get kallisto transcript counts from the 'PasillaTranscriptExpr' package
#' 
#' library(PasillaTranscriptExpr)
#' \donttest{
#' data_dir  <- system.file("extdata", package = "PasillaTranscriptExpr")
#' 
#' ## Load metadata
#' pasilla_metadata <- read.table(file.path(data_dir, "metadata.txt"), 
#' header = TRUE, as.is = TRUE)
#' 
#' ## Load counts
#' pasilla_counts <- read.table(file.path(data_dir, "counts.txt"), 
#' header = TRUE, as.is = TRUE)
#' 
#' ## Create a pasilla_samples data frame
#' pasilla_samples <- data.frame(sample_id = pasilla_metadata$SampleName, 
#'   group = pasilla_metadata$condition)
#' levels(pasilla_samples$group)
#' 
#' ## Create a dmDSdata object
#' d <- dmDSdata(counts = pasilla_counts, samples = pasilla_samples)
#' 
#' ## Use a subset of genes, which is defined in the following file
#' gene_id_subset <- readLines(file.path(data_dir, "gene_id_subset.txt"))
#' 
#' d <- d[names(d) %in% gene_id_subset, ]
#' 
#' # --------------------------------------------------------------------------
#' # Differential transcript usage analysis - simple two group comparison 
#' # --------------------------------------------------------------------------
#' 
#' ## Filtering
#' ## Check what is the minimal number of replicates per condition 
#' table(samples(d)$group)
#' 
#' d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3,
#'   min_gene_expr = 10, min_feature_expr = 10)
#' 
#' plotData(d)
#' 
#' ## Create the design matrix
#' design_full <- model.matrix(~ group, data = samples(d))
#' 
#' ## To make the analysis reproducible
#' set.seed(123)
#' ## Calculate precision
#' d <- dmPrecision(d, design = design_full)
#' 
#' plotPrecision(d)
#' 
#' head(mean_expression(d))
#' common_precision(d)
#' head(genewise_precision(d))
#' 
#' ## Fit full model proportions
#' d <- dmFit(d, design = design_full)
#' 
#' ## Get fitted proportions
#' head(proportions(d))
#' ## Get the DM regression coefficients (gene-level) 
#' head(coefficients(d))
#' ## Get the BB regression coefficients (feature-level) 
#' head(coefficients(d), level = "feature")
#' }
#' @author Malgorzata Nowicka
#' @seealso \code{\link{plotProportions}} \code{\link[edgeR]{glmFit}}
#' @references McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression 
#'   analysis of multifactor RNA-Seq experiments with respect to biological 
#'   variation. Nucleic Acids Research 40, 4288-4297.
#' @rdname dmFit
#' @importFrom limma nonEstimable
#' @export
setMethod("dmFit", "dmDSprecision", function(x, design, 
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
  
  if(!identical(x@design_precision, design))
    message(paste0("! The 'design' here is not identical as the 
      'design' used for precision estimation !\n"))
  
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
    precision = x@genewise_precision,
    one_way = one_way,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol,
    verbose = verbose, BPPARAM = BPPARAM)
  
  # Calculate the Beta-Binomial likelihoods for each feature
  if(bb_model){
    
    fit_bb <- bbDS_fit(counts = x@counts, fit = fit[["fit"]], design = design, 
      precision = x@genewise_precision,
      one_way = one_way,
      verbose = verbose, BPPARAM = BPPARAM)
    
    return(new("dmDSfit", design_fit_full = design, 
      fit_full = fit[["fit"]], 
      lik_full = fit[["lik"]], coef_full = fit[["coef"]],
      lik_full_bb = fit_bb[["lik"]], coef_full_bb = fit_bb[["coef"]],
      mean_expression = x@mean_expression, 
      common_precision = x@common_precision, 
      genewise_precision = x@genewise_precision, 
      design_precision = x@design_precision,
      counts = x@counts, samples = x@samples))
    
  }else{
    
    return(new("dmDSfit", design_fit_full = design, 
      fit_full = fit[["fit"]], 
      lik_full = fit[["lik"]], coef_full = fit[["coef"]],
      mean_expression = x@mean_expression, 
      common_precision = x@common_precision, 
      genewise_precision = x@genewise_precision, 
      design_precision = x@design_precision,
      counts = x@counts, samples = x@samples))
    
  }
  
  
})


################################################################################
### plotProportions
################################################################################

#' Plot feature proportions
#' 
#' This plot is available only for a group design, i.e., a design that is 
#' equivalent to multiple group fitting.
#' 
#' @return For a given gene, plot the observed and estimated with 
#'   Dirichlet-multinomial model feature proportions in each group. Estimated
#'   group proportions are marked with diamond shapes.
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
#' @param plot_fit Logical. Whether to plot the proportions estimated by the 
#'   full model.
#' @param plot_main Logical. Whether to plot a title with the information about 
#'   the Dirichlet-multinomial estimates.
#' @param group_colors Character vector with colors for each group defined by 
#'   \code{group_variable}.
#' @param feature_colors Character vector with colors for each feature of gene
#'   defined by \code{gene_id}.
#'   
#' @examples 
#' # --------------------------------------------------------------------------
#' # Create dmDSdata object 
#' # --------------------------------------------------------------------------
#' ## Get kallisto transcript counts from the 'PasillaTranscriptExpr' package
#' 
#' library(PasillaTranscriptExpr)
#' \donttest{
#' data_dir  <- system.file("extdata", package = "PasillaTranscriptExpr")
#' 
#' ## Load metadata
#' pasilla_metadata <- read.table(file.path(data_dir, "metadata.txt"), 
#' header = TRUE, as.is = TRUE)
#' 
#' ## Load counts
#' pasilla_counts <- read.table(file.path(data_dir, "counts.txt"), 
#' header = TRUE, as.is = TRUE)
#' 
#' ## Create a pasilla_samples data frame
#' pasilla_samples <- data.frame(sample_id = pasilla_metadata$SampleName, 
#'   group = pasilla_metadata$condition)
#' levels(pasilla_samples$group)
#' 
#' ## Create a dmDSdata object
#' d <- dmDSdata(counts = pasilla_counts, samples = pasilla_samples)
#' 
#' ## Use a subset of genes, which is defined in the following file
#' gene_id_subset <- readLines(file.path(data_dir, "gene_id_subset.txt"))
#' 
#' d <- d[names(d) %in% gene_id_subset, ]
#' 
#' # --------------------------------------------------------------------------
#' # Differential transcript usage analysis - simple two group comparison 
#' # --------------------------------------------------------------------------
#' 
#' ## Filtering
#' ## Check what is the minimal number of replicates per condition 
#' table(samples(d)$group)
#' 
#' d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3,
#'   min_gene_expr = 10, min_feature_expr = 10)
#' 
#' plotData(d)
#' 
#' ## Create the design matrix
#' design_full <- model.matrix(~ group, data = samples(d))
#' 
#' ## To make the analysis reproducible
#' set.seed(123)
#' ## Calculate precision
#' d <- dmPrecision(d, design = design_full)
#' 
#' plotPrecision(d)
#' 
#' head(mean_expression(d))
#' common_precision(d)
#' head(genewise_precision(d))
#' 
#' ## Fit full model proportions
#' d <- dmFit(d, design = design_full)
#' 
#' ## Get fitted proportions
#' head(proportions(d))
#' ## Get the DM regression coefficients (gene-level) 
#' head(coefficients(d))
#' ## Get the BB regression coefficients (feature-level) 
#' head(coefficients(d), level = "feature")
#' 
#' ## Fit null model proportions and perform the LR test to detect DTU
#' d <- dmTest(d, coef = "groupKD")
#' 
#' ## Plot the gene-level p-values
#' plotPValues(d)
#' 
#' ## Get the gene-level results
#' head(results(d))
#' 
#' ## Plot feature proportions for a top DTU gene
#' res <- results(d)
#' res <- res[order(res$pvalue, decreasing = FALSE), ]
#' 
#' top_gene_id <- res$gene_id[1]
#' 
#' plotProportions(d, gene_id = top_gene_id, group_variable = "group")
#' 
#' plotProportions(d, gene_id = top_gene_id, group_variable = "group", 
#'   plot_type = "lineplot")
#' 
#' plotProportions(d, gene_id = top_gene_id, group_variable = "group", 
#'   plot_type = "ribbonplot")
#' }
#' @author Malgorzata Nowicka
#' @seealso  \code{\link{plotData}}, \code{\link{plotPrecision}}, 
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
  
  if(is.factor(group))
    group <- factor(group)
  else
    group <- factor(group, levels = group)
  
  counts_gene <- x@counts[[gene_id]]
  
  if(!is.null(group_colors) && 
      plot_type %in% c("barplot", "boxplot1", "lineplot"))
    stopifnot(length(group_colors) == nlevels(group))
  if(!is.null(feature_colors) && 
      plot_type %in% c("boxplot2", "ribbonplot"))
    stopifnot(length(feature_colors) == nrow(counts_gene))
  
  if(nrow(counts_gene) <= 1)
    stop("!Gene has to have at least 2 features! \n")
  
  main <- NULL
  
  if(plot_main){
    
    mean_expression_gene <- mean(colSums(counts_gene), na.rm = TRUE)
    
    main <- paste0(gene_id, "\n Mean expression = ", 
      round(mean_expression_gene))
    
    precision_gene <- x@genewise_precision[gene_id]
    
    main <- paste0(main, ", Precision = ", round(precision_gene, 2))
    
  }
  
  prop_full <- NULL
  
  if(plot_fit){
    # Check if the design is equivalent to a one way layout
    groups <- edgeR::designAsFactor(x@design_fit_full)
    
    if(nlevels(groups) == ncol(x@design_fit_full) || 
        group_variable == "sample_id"){
      
      prop_full <- x@fit_full[[gene_id]][, !duplicated(group), drop = FALSE]
      colnames(prop_full) <- levels(group)
      
    }else{
      message("Fitted values are not plotted because the design does not 
        correspond to a group comparison defined by 'group_variable'!")
    }
    
  }
  
  # Order samples by group
  o <- order(group) 
  group <- group[o]
  counts_gene <- counts_gene[, o, drop = FALSE]
  
  ggp <- dm_plotProportions(counts = counts_gene, group = group, 
    prop_full = prop_full, main = main, plot_type = plot_type, 
    order = order, group_colors = group_colors, feature_colors = feature_colors)
  
  return(ggp)  
  
})


























