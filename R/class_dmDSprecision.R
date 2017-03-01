#' @include class_dmDSdata.R
NULL

################################################################################
### dmDSprecision class
################################################################################

#' dmDSprecision object
#' 
#' dmDSprecision extends the \code{\linkS4class{dmDSdata}} by adding the 
#' precision estimates of the Dirichlet-multinomial distribution used to model 
#' the feature (e.g., transcript, exon, exonic bin) counts for each gene in the 
#' differential usage analysis. Result of calling the \code{\link{dmPrecision}}
#' function.
#' 
#' @details Normally, in the differential analysis based on RNA-seq data, such 
#'   as, for example, differential gene expression, dispersion (of 
#'   negative-binomial model) is estimated. Here, we estimate precision of the 
#'   Dirichlet-multinomial model as it is more convenient computationally. To 
#'   obtain dispersion estimates, one can use a formula: dispersion = 1 / (1 + 
#'   precision).
#'   
#' @return
#' 
#' \itemize{ \item \code{mean_expression(x)}: Get a data frame with mean gene
#' expression. \item \code{common_precision(x), common_precision(x) <- value}:
#' Get or set common precision. \code{value} must be numeric of length 1. \item
#' \code{genewise_precision(x), genewise_precision(x) <- value}: Get a data
#' frame with gene-wise precision or set new gene-wise precision. \code{value}
#' must be a data frame with "gene_id" and "genewise_precision" columns. }
#' 
#' @param x,object dmDSprecision object.
#' @param value Values that replace current attributes.
#' @param ... Other parameters that can be defined by methods using this 
#'   generic.
#'   
#' @slot mean_expression Numeric vector of mean gene expression.
#' @slot common_precision Numeric value of estimated common precision.
#' @slot genewise_precision Numeric vector of estimated gene-wise precisions.
#' @slot design_precision Numeric matrix of the design used to estimate 
#'   precision.
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
#' }
#' @author Malgorzata Nowicka
#' @seealso \code{\linkS4class{dmDSdata}}, \code{\linkS4class{dmDSfit}}, 
#'   \code{\linkS4class{dmDStest}}
setClass("dmDSprecision", 
  contains = "dmDSdata",
  representation(mean_expression = "numeric", 
    common_precision = "numeric",
    genewise_precision = "numeric",
    design_precision = "matrix"))


# -----------------------------------------------------------------------------

setValidity("dmDSprecision", function(object){
  ## Has to return TRUE when valid object!
  out <- TRUE
  
  if(length(object@mean_expression) > 0){
    if(length(object@mean_expression) == length(object@counts)){
      if(all(names(object@mean_expression) == names(object@counts)))
        out <- TRUE
      else
        return(paste0("Different names of 'counts' and 'mean_expression'"))
    }
    else 
      return(paste0("Unequal length of 'counts' and 'mean_expression'"))
  }
  
  if(length(object@genewise_precision) > 0){
    if(length(object@genewise_precision) == length(object@counts)){
      if(all(names(object@genewise_precision) == names(object@counts)))
        out <- TRUE
      else
        return(paste0("Different names of 'counts' and 'genewise_precision'"))
    }
    else 
      return(paste0("Unequal length of 'counts' and 'genewise_precision'"))
  }
  
  if(length(object@common_precision) > 0){
    if(length(object@common_precision) == 1)
      out <- TRUE
    else
      return(paste0("'common_precision' must be a vector of length 1'"))
  }
  
  if(nrow(object@design_precision) == ncol(object@counts)){
    out <- TRUE
  }else{
    return(paste0("Number of rows in the design matrix must be equal 
          to the number of columns in counts"))
  }
  
  return(out)
  
})

################################################################################
### accessing methods
################################################################################

#' @rdname dmDSprecision-class
#' @param type Character indicating which design matrix should be returned.
#'   Possible values \code{"precision"}, \code{"full_model"} or
#'   \code{"null_model"}.
#' @export
setMethod("design", "dmDSprecision", function(object, type = "precision"){
  
  stopifnot(type %in% c("precision", "full_model", "null_model"))
  
  if(type == "precision")
    object@design_precision
  else
    NULL
  
})


#' @rdname dmDSprecision-class
#' @export
setGeneric("mean_expression", function(x, ...) 
  standardGeneric("mean_expression"))

#' @rdname dmDSprecision-class
#' @export
setMethod("mean_expression", "dmDSprecision", function(x){
  
  data.frame(gene_id = names(x@mean_expression), 
    mean_expression = x@mean_expression, 
    stringsAsFactors = FALSE, row.names = NULL)
  
})


#' @rdname dmDSprecision-class
#' @export
setGeneric("common_precision", function(x, ...) 
  standardGeneric("common_precision"))

#' @rdname dmDSprecision-class
#' @export
setMethod("common_precision", "dmDSprecision", function(x) 
  x@common_precision )


#' @rdname dmDSprecision-class
#' @export
setGeneric("common_precision<-", function(x, value) 
  standardGeneric("common_precision<-"))

#' @rdname dmDSprecision-class
#' @export
setMethod("common_precision<-", "dmDSprecision", function(x, value){
  ### value must be a numeric of length 1

  names(value) <- NULL
  
  return(new("dmDSprecision", mean_expression = x@mean_expression, 
    common_precision = value, genewise_precision = x@genewise_precision, 
    design_precision = x@design_precision,
    counts = x@counts, samples = x@samples))
  
})

#' @rdname dmDSprecision-class
#' @export
setGeneric("genewise_precision", function(x, ...) 
  standardGeneric("genewise_precision"))

#' @rdname dmDSprecision-class
#' @export
setMethod("genewise_precision", "dmDSprecision", function(x){
  
  data.frame(gene_id = names(x@genewise_precision), 
    genewise_precision = x@genewise_precision, stringsAsFactors = FALSE, 
    row.names = NULL)
  
})


#' @rdname dmDSprecision-class
#' @export
setGeneric("genewise_precision<-", function(x, value) 
  standardGeneric("genewise_precision<-"))


#' @rdname dmDSprecision-class
#' @export
setMethod("genewise_precision<-", "dmDSprecision", function(x, value){
  # value must be a data frame with gene_id and genewise_precision
  
  stopifnot(all(c("gene_id", "genewise_precision") %in% colnames(value)))
  stopifnot(all(names(x@counts) %in% value[,"gene_id"]))
  order <- match(names(x@counts), value[,"gene_id"])
  
  return(new("dmDSprecision", mean_expression = x@mean_expression, 
    common_precision = x@common_precision, 
    genewise_precision = value[order, "genewise_precision"],
    design_precision = x@design_precision,  
    counts = x@counts, samples = x@samples))
  
})



# -----------------------------------------------------------------------------


setMethod("show", "dmDSprecision", function(object){
  
  callNextMethod(object)
  
  cat("  design()\n")
  cat("  mean_expression(), common_precision(), genewise_precision()\n")

})


################################################################################
### dmPrecision
################################################################################

#' Estimate the precision parameter in the Dirichlet-multinomial model
#' 
#' Maximum likelihood estimates of the precision parameter in the 
#' Dirichlet-multinomial model used for the differential exon/transcript usage
#' or QTL analysis.
#' 
#' @details Normally, in the differential analysis based on RNA-seq data, such 
#'   as, for example, differential gene expression, dispersion (of 
#'   negative-binomial model) is estimated. Here, we estimate precision of the 
#'   Dirichlet-multinomial model as it is more convenient computationally. To 
#'   obtain dispersion estimates, one can use a formula: dispersion = 1 / (1 + 
#'   precision).
#'   
#' @param x \code{\linkS4class{dmDSdata}} or \code{\linkS4class{dmSQTLdata}} 
#'   object.
#' @param ... Other parameters that can be defined by methods using this 
#'   generic.
#' @export
setGeneric("dmPrecision", function(x, ...) standardGeneric("dmPrecision"))


# -----------------------------------------------------------------------------


#' @details Parameters that are used in the precision (dispersion = 1 / (1 + 
#'   precision)) estimation start with prefix \code{prec_}. Those that are used 
#'   for the proportion estimation in each group when the shortcut fitting 
#'   \code{one_way = TRUE} can be used start with \code{prop_}, and those that 
#'   are used in the regression framework start with \code{coef_}.
#'   
#'   There are two optimization methods implemented within dmPrecision: 
#'   \code{"optimize"} for the common precision and \code{"grid"} for the 
#'   gene-wise precision.
#'   
#'   Only part of the precision parameters in dmPrecision have an influence on
#'   a given optimization method. Here is a list of such active parameters:
#'   
#'   \code{"optimize"}:
#'   
#'   \itemize{ \item \code{prec_interval}: Passed as \code{interval}. \item 
#'   \code{prec_tol}: The accuracy defined as \code{tol}. }
#'   
#'   \code{"grid"}, which uses the grid approach from 
#'   \code{\link[edgeR]{estimateDisp}} in \code{\link{edgeR}}:
#'   
#'   \itemize{ \item \code{prec_init}, \code{prec_grid_length}, 
#'   \code{prec_grid_range}: Parameters used to construct the search grid 
#'   \code{prec_init * 2^seq(from = prec_grid_range[1]}, \code{to = 
#'   prec_grid_range[2]}, \code{length = prec_grid_length)}. \item 
#'   \code{prec_moderation}: Dipsersion shrinkage is available only with 
#'   \code{"grid"} method. \item \code{prec_prior_df}: Used only when precision
#'   shrinkage is activated. Moderated likelihood is equal to \code{loglik + 
#'   prec_prior_df * moderation}. Higher \code{prec_prior_df}, more shrinkage 
#'   toward common or trended precision is applied. \item \code{prec_span}: 
#'   Used only when precision moderation toward trend is activated. }
#'   
#' @param design Numeric matrix defining the model that should be used when 
#'   estimating precision. Normally this should be a full model design used 
#'   also in \code{\link{dmFit}}.
#' @param mean_expression Logical. Whether to estimate the mean expression of 
#'   genes.
#' @param common_precision Logical. Whether to estimate the common precision.
#' @param genewise_precision Logical. Whether to estimate the gene-wise 
#'   precision.
#' @param prec_adjust Logical. Whether to use the Cox-Reid adjusted or 
#'   non-adjusted profile likelihood.
#' @param one_way Logical. Should the shortcut fitting be used when the design 
#'   corresponds to multiple group comparison. This is a similar approach as in 
#'   \code{\link{edgeR}}. If \code{TRUE} (the default), then proportions are 
#'   fitted per group and regression coefficients are recalculated from those 
#'   fits.
#' @param prec_subset Value from 0 to 1 defining the percentage of genes used in
#'   common precision estimation. The default is 0.1, which uses 10% of 
#'   randomly selected genes to speed up the precision estimation process. Use 
#'   \code{set.seed} function to make the analysis reproducible. See Examples.
#' @param prec_interval Numeric vector of length 2 defining the interval of 
#'   possible values for the common precision.
#' @param prec_tol The desired accuracy when estimating common precision.
#' @param prec_init Initial precision. If \code{common_precision} is 
#'   \code{TRUE}, then \code{prec_init} is overwritten by common precision 
#'   estimate.
#' @param prec_grid_length Length of the search grid.
#' @param prec_grid_range Vector giving the limits of grid interval.
#' @param prec_moderation Precision moderation method. One can choose to shrink
#'   the precision estimates toward the common precision (\code{"common"}) or 
#'   toward the (precision versus mean expression) trend (\code{"trended"})
#' @param prec_prior_df Degree of moderation (shrinkage) in case when it can not
#'   be calculated automaticaly (number of genes on the upper boundary of grid 
#'   is smaller than 10). By default it is equal to 0.
#' @param prec_span Value from 0 to 1 defining the percentage of genes used in 
#'   smoothing sliding window when calculating the precision versus mean 
#'   expression trend.
#' @param prop_mode Optimization method used to estimate proportions. Possible 
#'   value \code{"constrOptim"}.
#' @param prop_tol The desired accuracy when estimating proportions.
#' @param coef_mode Optimization method used to estimate regression 
#'   coefficients. Possible value \code{"optim"}.
#' @param coef_tol The desired accuracy when estimating regression coefficients.
#' @param verbose Numeric. Definie the level of progress messages displayed. 0 -
#'   no messages, 1 - main messages, 2 - message for every gene fitting.
#' @param BPPARAM Parallelization method used by 
#'   \code{\link[BiocParallel]{bplapply}}.
#'   
#' @return Returns a \code{\linkS4class{dmDSprecision}} or 
#'   \code{\linkS4class{dmSQTLprecision}} object.
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
#' }
#' @seealso \code{\link{plotPrecision}} \code{\link[edgeR]{estimateDisp}}
#' @author Malgorzata Nowicka
#' @references McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression
#' analysis of multifactor RNA-Seq experiments with respect to biological
#' variation. Nucleic Acids Research 40, 4288-4297.
#' @rdname dmPrecision
#' @importFrom limma nonEstimable
#' @export
setMethod("dmPrecision", "dmDSdata", function(x, design, 
  mean_expression = TRUE, common_precision = TRUE, genewise_precision = TRUE,
  prec_adjust = TRUE, prec_subset = 0.1, 
  prec_interval = c(0, 1e+3), prec_tol = 1e+01, 
  prec_init = 100, prec_grid_length = 21, prec_grid_range = c(-10, 10), 
  prec_moderation = "trended", prec_prior_df = 0, prec_span = 0.1, 
  one_way = TRUE, 
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
  
  # Check other parameters
  stopifnot(is.logical(mean_expression))
  stopifnot(is.logical(common_precision))
  stopifnot(is.logical(genewise_precision))
  stopifnot(is.logical(prec_adjust))
  stopifnot(length(prec_subset) == 1)
  stopifnot(is.numeric(prec_subset) && prec_subset > 0 && prec_subset <= 1)
  stopifnot(length(prec_interval) == 2)
  stopifnot(prec_interval[1] < prec_interval[2])
  stopifnot(length(prec_tol) == 1)
  stopifnot(is.numeric(prec_tol) && prec_tol > 0)
  stopifnot(length(prec_init) == 1)
  stopifnot(is.numeric(prec_init))
  stopifnot(prec_grid_length > 2)
  stopifnot(length(prec_grid_range) == 2)
  stopifnot(prec_grid_range[1] < prec_grid_range[2])
  stopifnot(length(prec_moderation) == 1)
  stopifnot(prec_moderation %in% c("none", "common", "trended"))
  stopifnot(length(prec_prior_df) == 1)
  stopifnot(is.numeric(prec_prior_df) && prec_prior_df >= 0)
  stopifnot(length(prec_span) == 1)
  stopifnot(is.numeric(prec_span) && prec_span > 0 && prec_span < 1)
  
  stopifnot(is.logical(one_way))
  stopifnot(length(prop_mode) == 1)
  stopifnot(prop_mode %in% c("constrOptim"))
  stopifnot(length(prop_tol) == 1)
  stopifnot(is.numeric(prop_tol) && prop_tol > 0)
  
  stopifnot(length(coef_mode) == 1)
  stopifnot(coef_mode %in% c("optim"))
  stopifnot(length(coef_tol) == 1)
  stopifnot(is.numeric(coef_tol) && coef_tol > 0)
  
  stopifnot(verbose %in% 0:2)
  
  if(mean_expression || (genewise_precision && 
      prec_moderation == "trended")){
    mean_expression <- dm_estimateMeanExpression(counts = x@counts, 
      verbose = verbose)
  }else{
    mean_expression <- numeric()
  }
  
  if(common_precision){
    
    if(prec_subset < 1){
      
      message(paste0("! Using a subset of ", prec_subset, 
        " genes to estimate common precision !\n"))
      
      genes2keep <- sample(1:length(x@counts), 
        max(round(prec_subset * length(x@counts)), 1), replace = FALSE)
      
    }else{
      genes2keep <- 1:length(x@counts)
    }
    
    common_precision <- dmDS_estimateCommonPrecision(
      counts = x@counts[genes2keep, ], 
      design = design, prec_adjust = prec_adjust, 
      prec_interval = prec_interval, prec_tol = prec_tol,
      one_way = one_way,
      prop_mode = prop_mode, prop_tol = prop_tol, 
      coef_mode = coef_mode, coef_tol = coef_tol,
      verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    common_precision <- numeric()
  }
  
  if(genewise_precision){
    
    if(length(common_precision) == 1){
      message("! Using common_precision = ", round(common_precision, 4), 
        " as prec_init !\n")
      prec_init <- common_precision
    }
    
    genewise_precision <- dmDS_estimateTagwisePrecision(counts = x@counts, 
      design = design, mean_expression = mean_expression, 
      prec_adjust = prec_adjust, prec_init = prec_init, 
      prec_grid_length = prec_grid_length, prec_grid_range = prec_grid_range, 
      prec_moderation = prec_moderation, 
      prec_prior_df = prec_prior_df, prec_span = prec_span, 
      one_way = one_way,
      prop_mode = prop_mode, prop_tol = prop_tol, 
      coef_mode = coef_mode, coef_tol = coef_tol,
      verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    genewise_precision <- numeric()
  }
  
  return(new("dmDSprecision", mean_expression = mean_expression, 
    common_precision = common_precision, 
    genewise_precision = genewise_precision, 
    design_precision = design,
    counts = x@counts, samples = x@samples))
  
})


################################################################################
### plotPrecision
################################################################################

#' Precision versus mean expression plot
#' 
#' @return Normally in the differential analysis based on RNA-seq data, such 
#'   plot has dispersion parameter plotted on the y-axis. Here, the y-axis 
#'   represents precision since in the Dirichlet-multinomial model this is the 
#'   parameter that is directly estimated. It is important to keep in mind that
#'   the precision parameter (gamma0) is inverse proportional to dispersion 
#'   (theta): theta = 1 / (1 + gamma0). In RNA-seq data, we can typically
#'   observe a trend where the dispersion decreases (here, precision increases)
#'   for genes with higher mean expression.
#'   
#' @param x \code{\linkS4class{dmDSprecision}} or 
#'   \code{\linkS4class{dmSQTLprecision}} object.
#' @param ... Other parameters that can be defined by methods using this 
#'   generic.
#' @export
setGeneric("plotPrecision", function(x, ...) standardGeneric("plotPrecision"))


# -----------------------------------------------------------------------------


#' @inheritParams plotData
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
#' }
#' @author Malgorzata Nowicka
#' @seealso \code{\link{plotData}}, \code{\link{plotProportions}},
#'   \code{\link{plotPValues}}
#'   
#' @rdname plotPrecision
#' @export
setMethod("plotPrecision", "dmDSprecision", function(x){
  
  if(!length(x@genewise_precision) == length(x@counts))
    stop("Genewise precision must be estimated for each gene!")
  if(!length(x@genewise_precision) == length(x@mean_expression))
    stop("Mean expression must be estimated for each gene!")
  
  if(length(x@common_precision) == 0){
    common_precision <- NULL
  }else{
    common_precision <- x@common_precision
  }
  
  ggp <- dm_plotPrecision(genewise_precision = x@genewise_precision, 
    mean_expression = x@mean_expression, nr_features = elementNROWS(x@counts), 
    common_precision = common_precision)
  
  return(ggp)  
  
})








