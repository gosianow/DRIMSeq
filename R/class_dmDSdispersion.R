#' @include class_dmDSdata.R
NULL

################################################################################
### dmDSdispersion class
################################################################################

#' dmDSdispersion object
#' 
#' dmDSdispersion extends the \code{\linkS4class{dmDSdata}} by adding the 
#' precision estimates of the Dirichlet-multinomial distribution used to model 
#' the feature (e.g., transcript, exon, exonic bin) counts for each gene in the 
#' differential usage analysis. Result of calling the \code{\link{dmDispersion}}
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
#' expression. \item \code{common_dispersion(x), common_dispersion(x) <- value}:
#' Get or set common dispersion. \code{value} must be numeric of length 1. \item
#' \code{genewise_dispersion(x), genewise_dispersion(x) <- value}: Get a data
#' frame with gene-wise dispersion or set new gene-wise dispersion. \code{value}
#' must be a data frame with "gene_id" and "genewise_dispersion" columns. }
#' 
#' @param x dmDSdispersion object.
#' @param value Values that replace current attributes.
#' @param ... Other parameters that can be defined by methods using this 
#'   generic.
#'   
#' @slot mean_expression Numeric vector of mean gene expression.
#' @slot common_dispersion Numeric value of estimated common dispersion.
#' @slot genewise_dispersion Numeric vector of estimated gene-wise dispersions.
#' @slot design_dispersion Numeric matrix of the desing used to estimate 
#'   dispersion.
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
#' metadata <- read.table(file.path(data_dir, "metadata.txt"), header = TRUE, 
#'   as.is = TRUE)
#' 
#' ## Load counts
#' counts <- read.table(file.path(data_dir, "counts.txt"), header = TRUE, 
#'   as.is = TRUE)
#' 
#' ## Create a samples data frame
#' samples <- data.frame(sample_id = metadata$SampleName, 
#'   group = metadata$condition)
#' levels(samples$group)
#' 
#' ## Create a dmDSdata object
#' d <- dmDSdata(counts = counts, samples = samples)
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
#' design <- model.matrix(~ group, data = samples(d))
#' 
#' ## To make the analysis reproducible
#' set.seed(123)
#' ## Calculate dispersion
#' d <- dmDispersion(d, design = design)
#' 
#' plotDispersion(d)
#' 
#' head(mean_expression(d))
#' common_dispersion(d)
#' head(genewise_dispersion(d))
#' }
#' @author Malgorzata Nowicka
#' @seealso \code{\linkS4class{dmDSdata}}, \code{\linkS4class{dmDSfit}}, 
#'   \code{\linkS4class{dmDStest}}
setClass("dmDSdispersion", 
  contains = "dmDSdata",
  representation(mean_expression = "numeric", 
    common_dispersion = "numeric",
    genewise_dispersion = "numeric",
    design_dispersion = "matrix"))


# -----------------------------------------------------------------------------

setValidity("dmDSdispersion", function(object){
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
  
  if(length(object@genewise_dispersion) > 0){
    if(length(object@genewise_dispersion) == length(object@counts)){
      if(all(names(object@genewise_dispersion) == names(object@counts)))
        out <- TRUE
      else
        return(paste0("Different names of 'counts' and 'genewise_dispersion'"))
    }
    else 
      return(paste0("Unequal length of 'counts' and 'genewise_dispersion'"))
  }
  
  if(length(object@common_dispersion) > 0){
    if(length(object@common_dispersion) == 1)
      out <- TRUE
    else
      return(paste0("'common_dispersion' must be a vector of length 1'"))
  }
  
  if(nrow(object@design_dispersion) == ncol(object@counts)){
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


#' @rdname dmDSdispersion-class
#' @export
setGeneric("mean_expression", function(x, ...) 
  standardGeneric("mean_expression"))

#' @rdname dmDSdispersion-class
#' @export
setMethod("mean_expression", "dmDSdispersion", function(x){
  
  data.frame(gene_id = names(x@mean_expression), 
    mean_expression = x@mean_expression, 
    stringsAsFactors = FALSE, row.names = NULL)
  
})


#' @rdname dmDSdispersion-class
#' @export
setGeneric("common_dispersion", function(x, ...) 
  standardGeneric("common_dispersion"))

#' @rdname dmDSdispersion-class
#' @export
setMethod("common_dispersion", "dmDSdispersion", function(x) 
  x@common_dispersion )


#' @rdname dmDSdispersion-class
#' @export
setGeneric("common_dispersion<-", function(x, value) 
  standardGeneric("common_dispersion<-"))

#' @rdname dmDSdispersion-class
#' @export
setMethod("common_dispersion<-", "dmDSdispersion", function(x, value){
  ### value must be a numeric of length 1

  names(value) <- NULL
  
  return(new("dmDSdispersion", mean_expression = x@mean_expression, 
    common_dispersion = value, genewise_dispersion = x@genewise_dispersion, 
    design_dispersion = x@design_dispersion,
    counts = x@counts, samples = x@samples))
  
})

#' @rdname dmDSdispersion-class
#' @export
setGeneric("genewise_dispersion", function(x, ...) 
  standardGeneric("genewise_dispersion"))

#' @rdname dmDSdispersion-class
#' @export
setMethod("genewise_dispersion", "dmDSdispersion", function(x){
  
  data.frame(gene_id = names(x@genewise_dispersion), 
    genewise_dispersion = x@genewise_dispersion, stringsAsFactors = FALSE, 
    row.names = NULL)
  
})


#' @rdname dmDSdispersion-class
#' @export
setGeneric("genewise_dispersion<-", function(x, value) 
  standardGeneric("genewise_dispersion<-"))


#' @rdname dmDSdispersion-class
#' @export
setMethod("genewise_dispersion<-", "dmDSdispersion", function(x, value){
  # value must be a data frame with gene_id and genewise_dispersion
  
  stopifnot(all(c("gene_id", "genewise_dispersion") %in% colnames(value)))
  stopifnot(all(names(x@counts) %in% value[,"gene_id"]))
  order <- match(names(x@counts), value[,"gene_id"])
  
  return(new("dmDSdispersion", mean_expression = x@mean_expression, 
    common_dispersion = x@common_dispersion, 
    genewise_dispersion = value[order, "genewise_dispersion"],
    design_dispersion = x@design_dispersion,  
    counts = x@counts, samples = x@samples))
  
})



# -----------------------------------------------------------------------------


setMethod("show", "dmDSdispersion", function(object){
  
  callNextMethod(object)
  
  cat("  mean_expression(), common_dispersion(), genewise_dispersion()\n")
  
})


################################################################################
### dmDispersion
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
setGeneric("dmDispersion", function(x, ...) standardGeneric("dmDispersion"))


# -----------------------------------------------------------------------------


#' @details Parameters that are used in the precision (dispersion = 1 / (1 + 
#'   precision)) estimation start with prefix \code{disp_}. Those that are used 
#'   for the proportion estimation in each group when the shortcut fitting 
#'   \code{one_way = TRUE} can be used start with \code{prop_}, and those that 
#'   are used in the regression framework start with \code{coef_}.
#'   
#'   There are two optimization methods implemented within dmDispersion: 
#'   \code{"optimize"} for the common dispersion and \code{"grid"} for the 
#'   gene-wise dispersion.
#'   
#'   Only part of the dispersion parameters in dmDispersion have an influence on
#'   a given optimization method. Here is a list of such active parameters:
#'   
#'   \code{"optimize"}:
#'   
#'   \itemize{ \item \code{disp_interval}: Passed as \code{interval}. \item 
#'   \code{disp_tol}: The accuracy defined as \code{tol}. }
#'   
#'   \code{"grid"}, which uses the grid approach from 
#'   \code{\link[edgeR]{estimateDisp}} in \code{\link{edgeR}}:
#'   
#'   \itemize{ \item \code{disp_init}, \code{disp_grid_length}, 
#'   \code{disp_grid_range}: Parameters used to construct the search grid 
#'   \code{disp_init * 2^seq(from = disp_grid_range[1]}, \code{to = 
#'   disp_grid_range[2]}, \code{length = disp_grid_length)}. \item 
#'   \code{disp_moderation}: Dipsersion shrinkage is available only with 
#'   \code{"grid"} method. \item \code{disp_prior_df}: Used only when dispersion
#'   shrinkage is activated. Moderated likelihood is equal to \code{loglik + 
#'   disp_prior_df * moderation}. Higher \code{disp_prior_df}, more shrinkage 
#'   toward common or trended dispersion is applied. \item \code{disp_span}: 
#'   Used only when dispersion moderation toward trend is activated. }
#'   
#' @param design Numeric matrix definig the model that should be used when 
#'   estimating dispersion. Normally this should be a full model design used 
#'   also in \code{\link{dmFit}}.
#' @param mean_expression Logical. Whether to estimate the mean expression of 
#'   genes.
#' @param common_dispersion Logical. Whether to estimate the common dispersion.
#' @param genewise_dispersion Logical. Whether to estimate the gene-wise 
#'   dispersion.
#' @param disp_adjust Logical. Whether to use the Cox-Reid adjusted or 
#'   non-adjusted profile likelihood.
#' @param one_way Logical. Should the shortcut fitting be used when the design 
#'   corresponds to multiple group comparison. This is a similar approach as in 
#'   \code{\link{edgeR}}. If \code{TRUE} (the default), then proportions are 
#'   fitted per group and regression coefficients are recalculated from those 
#'   fits.
#' @param disp_subset Value from 0 to 1 defining the percentage of genes used in
#'   common dispersion estimation. The default is 0.1, which uses 10% of 
#'   randomly selected genes to speed up the dispersion estimation process. Use 
#'   \code{set.seed} function to make the analysis reproducible. See Examples.
#' @param disp_interval Numeric vector of length 2 defining the interval of 
#'   possible values for the common dispersion.
#' @param disp_tol The desired accuracy when estimating common dispersion.
#' @param disp_init Initial dispersion. If \code{common_dispersion} is 
#'   \code{TRUE}, then \code{disp_init} is overwritten by common dispersion 
#'   estimate.
#' @param disp_grid_length Length of the search grid.
#' @param disp_grid_range Vector giving the limits of grid interval.
#' @param disp_moderation Dispersion moderation method. One can choose to shrink
#'   the dispersion estimates toward the common dispersion (\code{"common"}) or 
#'   toward the (dispersion versus mean expression) trend (\code{"trended"})
#' @param disp_prior_df Degree of moderation (shrinkage) in case when it can not
#'   be calculated automaticaly (number of genes on the upper boundary of grid 
#'   is smaller than 10). By default it is equal to 0.
#' @param disp_span Value from 0 to 1 defining the percentage of genes used in 
#'   smoothing sliding window when calculating the dispersion versus mean 
#'   expression trend.
#' @param prop_mode Optimization method used to estimate proportions. Possible 
#'   value \code{"constrOptim"}.
#' @param prop_tol The desired accuracy when estimating proportions.
#' @param coef_mode Optimization method used to estimate regression 
#'   coefficients. Possible values \code{"optim"} (the default), \code{"nlminb"}
#'   or \code{"Rcgmin"}.
#' @param coef_tol The desired accuracy when estimating regression coefficients.
#' @param verbose Numeric. Definie the level of progress messages displayed. 0 -
#'   no messages, 1 - main messages, 2 - message for every gene fitting.
#' @param BPPARAM Parallelization method used by 
#'   \code{\link[BiocParallel]{bplapply}}.
#'   
#' @return Returns a \code{\linkS4class{dmDSdispersion}} or 
#'   \code{\linkS4class{dmSQTLdispersion}} object.
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
#' metadata <- read.table(file.path(data_dir, "metadata.txt"), header = TRUE, 
#'   as.is = TRUE)
#' 
#' ## Load counts
#' counts <- read.table(file.path(data_dir, "counts.txt"), header = TRUE, 
#'   as.is = TRUE)
#' 
#' ## Create a samples data frame
#' samples <- data.frame(sample_id = metadata$SampleName, 
#'   group = metadata$condition)
#' levels(samples$group)
#' 
#' ## Create a dmDSdata object
#' d <- dmDSdata(counts = counts, samples = samples)
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
#' design <- model.matrix(~ group, data = samples(d))
#' 
#' ## To make the analysis reproducible
#' set.seed(123)
#' ## Calculate dispersion
#' d <- dmDispersion(d, design = design)
#' 
#' plotDispersion(d)
#' 
#' head(mean_expression(d))
#' common_dispersion(d)
#' head(genewise_dispersion(d))
#' }
#' @seealso \code{\link{plotDispersion}} \code{\link[edgeR]{estimateDisp}}
#' @author Malgorzata Nowicka
#' @references McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression
#' analysis of multifactor RNA-Seq experiments with respect to biological
#' variation. Nucleic Acids Research 40, 4288-4297.
#' @rdname dmDispersion
#' @export
setMethod("dmDispersion", "dmDSdata", function(x, design, 
  mean_expression = TRUE, common_dispersion = TRUE, genewise_dispersion = TRUE,
  disp_adjust = TRUE, disp_subset = 0.1, 
  disp_interval = c(0, 1e+5), disp_tol = 1e+01, 
  disp_init = 100, disp_grid_length = 21, disp_grid_range = c(-10, 10), 
  disp_moderation = "trended", disp_prior_df = 0, disp_span = 0.1, 
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
  stopifnot(is.logical(common_dispersion))
  stopifnot(is.logical(genewise_dispersion))
  stopifnot(is.logical(disp_adjust))
  stopifnot(length(disp_subset) == 1)
  stopifnot(is.numeric(disp_subset) && disp_subset > 0 && disp_subset <= 1)
  stopifnot(length(disp_interval) == 2)
  stopifnot(disp_interval[1] < disp_interval[2])
  stopifnot(length(disp_tol) == 1)
  stopifnot(is.numeric(disp_tol) && disp_tol > 0)
  stopifnot(length(disp_init) == 1)
  stopifnot(is.numeric(disp_init))
  stopifnot(disp_grid_length > 2)
  stopifnot(length(disp_grid_range) == 2)
  stopifnot(disp_grid_range[1] < disp_grid_range[2])
  stopifnot(length(disp_moderation) == 1)
  stopifnot(disp_moderation %in% c("none", "common", "trended"))
  stopifnot(length(disp_prior_df) == 1)
  stopifnot(is.numeric(disp_prior_df) && disp_prior_df >= 0)
  stopifnot(length(disp_span) == 1)
  stopifnot(is.numeric(disp_span) && disp_span > 0 && disp_span < 1)
  
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
  
  if(mean_expression || (genewise_dispersion && 
      disp_moderation == "trended")){
    mean_expression <- dm_estimateMeanExpression(counts = x@counts, 
      verbose = verbose)
  }else{
    mean_expression <- numeric()
  }
  
  if(common_dispersion){
    
    if(disp_subset < 1){
      
      message(paste0("! Using a subset of ", disp_subset, 
        " genes to estimate common dispersion !\n"))
      
      genes2keep <- sample(1:length(x@counts), 
        max(round(disp_subset * length(x@counts)), 1), replace = FALSE)
      
    }else{
      genes2keep <- 1:length(x@counts)
    }
    
    common_dispersion <- dmDS_estimateCommonDispersion(
      counts = x@counts[genes2keep, ], 
      design = design, disp_adjust = disp_adjust, 
      disp_interval = disp_interval, disp_tol = disp_tol,
      one_way = one_way,
      prop_mode = prop_mode, prop_tol = prop_tol, 
      coef_mode = coef_mode, coef_tol = coef_tol,
      verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    common_dispersion <- numeric()
  }
  
  if(genewise_dispersion){
    
    if(length(common_dispersion) == 1){
      message("! Using common_dispersion = ", round(common_dispersion, 4), 
        " as disp_init !\n")
      disp_init <- common_dispersion
    }
    
    genewise_dispersion <- dmDS_estimateTagwiseDispersion(counts = x@counts, 
      design = design, mean_expression = mean_expression, 
      disp_adjust = disp_adjust, disp_init = disp_init, 
      disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, 
      disp_moderation = disp_moderation, 
      disp_prior_df = disp_prior_df, disp_span = disp_span, 
      one_way = one_way,
      prop_mode = prop_mode, prop_tol = prop_tol, 
      coef_mode = coef_mode, coef_tol = coef_tol,
      verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    genewise_dispersion <- numeric()
  }
  
  return(new("dmDSdispersion", mean_expression = mean_expression, 
    common_dispersion = common_dispersion, 
    genewise_dispersion = genewise_dispersion, 
    design_dispersion = design,
    counts = x@counts, samples = x@samples))
  
})


################################################################################
### plotDispersion
################################################################################

#' Precision versus mean expression plot
#' 
#' @return Normally in the differential analysis based on RNA-seq data, such 
#'   plot has dispersion parameter plotted on the y-axis. Here, the y-axis 
#'   represents precision since in the Dirichlet-multinomial model this is the 
#'   parameter that is directly estimated. It is important to keep in mind that
#'   the precision parameter (gamma0) is inverse proportional to dispersion 
#'   (theta): theta = 1 / (1 + gamma0). In RNA-seq data, we can typically
#'   observe a trend where the dispersion decreases (here, presicion increases)
#'   for genes with higher mean expression.
#'   
#' @param x \code{\linkS4class{dmDSdispersion}} or 
#'   \code{\linkS4class{dmSQTLdispersion}} object.
#' @param ... Other parameters that can be defined by methods using this 
#'   generic.
#' @export
setGeneric("plotDispersion", function(x, ...) standardGeneric("plotDispersion"))


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
#' metadata <- read.table(file.path(data_dir, "metadata.txt"), header = TRUE, 
#'   as.is = TRUE)
#' 
#' ## Load counts
#' counts <- read.table(file.path(data_dir, "counts.txt"), header = TRUE, 
#'   as.is = TRUE)
#' 
#' ## Create a samples data frame
#' samples <- data.frame(sample_id = metadata$SampleName, 
#'   group = metadata$condition)
#' levels(samples$group)
#' 
#' ## Create a dmDSdata object
#' d <- dmDSdata(counts = counts, samples = samples)
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
#' design <- model.matrix(~ group, data = samples(d))
#' 
#' ## To make the analysis reproducible
#' set.seed(123)
#' ## Calculate dispersion
#' d <- dmDispersion(d, design = design)
#' 
#' plotDispersion(d)
#' 
#' head(mean_expression(d))
#' common_dispersion(d)
#' head(genewise_dispersion(d))
#' }
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmDSdata}}, \code{\link{data_dmSQTLdata}},
#'   \code{\link{plotData}}, \code{\link{plotProportions}}, \code{\link{plotPValues}}
#'   
#' @rdname plotDispersion
#' @export
setMethod("plotDispersion", "dmDSdispersion", function(x){
  
  if(!length(x@genewise_dispersion) == length(x@counts))
    stop("Genewise dispersion must be estimated for each gene!")
  if(!length(x@genewise_dispersion) == length(x@mean_expression))
    stop("Mean expression must be estimated for each gene!")
  
  if(length(x@common_dispersion) == 0){
    common_dispersion <- NULL
  }else{
    common_dispersion <- x@common_dispersion
  }
  
  ggp <- dm_plotDispersion(genewise_dispersion = x@genewise_dispersion, 
    mean_expression = x@mean_expression, nr_features = elementNROWS(x@counts), 
    common_dispersion = common_dispersion)
  
  return(ggp)  
  
})








