#' @include class_dmDSdata.R
NULL

################################################################################
### dmDSdispersion class
################################################################################

#' dmDSdispersion object
#' 
#' dmDSdispersion extends the \code{\linkS4class{dmDSdata}} by adding the
#' dispersion estimates of Dirichlet-multinomial distribution used to model the
#' feature (e.g., transcript, exon, exonic bin) coiunts for each gene in the
#' differential splicing analysis. Result of \code{\link{dmDispersion}}.
#' 
#' @return
#' 
#' \itemize{ 
#' \item \code{mean_expression(x)}: Get a data frame with mean gene
#' expression. 
#' \item \code{common_dispersion(x), common_dispersion(x) <- value}:
#' Get or set common dispersion. \code{value} must be numeric of length 1. 
#' \item \code{genewise_dispersion(x), genewise_dispersion(x) <- value}: 
#' Get a data frame with gene-wise dispersion or set new gene-wise dispersion. 
#' \code{value} must be a data frame with "gene_id" and "genewise_dispersion"
#' columns. 
#' }
#' 
#' 
#' @param x dmDSdispersion object.
#' @param value Values that replace current attributes. 
#' @param ... Other parameters that can be defined by methods using this
#'   generic.
#'   
#' @slot mean_expression Numeric vector of mean gene expression.
#' @slot common_dispersion Numeric value of estimated common dispersion.
#' @slot genewise_dispersion Numeric vector of estimated gene-wise dispersions.
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
#' plotDispersion(d)
#' 
#' head(mean_expression(d))
#' common_dispersion(d)
#' head(genewise_dispersion(d))
#' 
#' }
#' 
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmDSdata}}, \code{\linkS4class{dmDSdata}},
#'   \code{\linkS4class{dmDSfit}}, \code{\linkS4class{dmDStest}}
setClass("dmDSdispersion", 
  contains = "dmDSdata",
  representation(mean_expression = "numeric", common_dispersion = "numeric",
    genewise_dispersion = "numeric"))


####################################

setValidity("dmDSdispersion", function(object){
  # has to return TRUE when valid object!
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
  ### value must be a data.frame with gene_id and genewise_dispersion columns
  stopifnot(all(c("gene_id", "genewise_dispersion") %in% colnames(value)))
  stopifnot(all(names(x@counts) %in% value[,"gene_id"]))
  order <- match(names(x@counts), value[,"gene_id"])
  
  return(new("dmDSdispersion", mean_expression = x@mean_expression, 
    common_dispersion = x@common_dispersion, 
    genewise_dispersion = value[order, "genewise_dispersion"], 
    counts = x@counts, samples = x@samples))
  
})



#################################

setMethod("show", "dmDSdispersion", function(object){
  
  callNextMethod(object)
  
  cat("  mean_expression(), common_dispersion(), genewise_dispersion()\n")
  
})


################################################################################
### dmDispersion
################################################################################

#' Estimate dispersions in Dirichlet-multinomial model
#' 
#' Maximum likelihood estimates of dispersion parameters in the
#' Dirichlet-multinomial model used in differential splicing or sQTL analysis.
#' 
#' @param x \code{\linkS4class{dmDSdata}} or \code{\linkS4class{dmSQTLdata}}
#'   object.
#' @param ... Other parameters that can be defined by methods using this
#'   generic.
#' @export
setGeneric("dmDispersion", function(x, ...) standardGeneric("dmDispersion"))


#################################

#' @details Parameters that are used in the dispersion estimation start with 
#'   prefix \code{disp_}, and those that are used for the proportion estimation 
#'   start with \code{prop_}.
#'   
#'   There are 4 optimization methods implemented within dmDispersion 
#'   (\code{"optimize"}, \code{"optim"}, \code{"constrOptim"} and \code{"grid"})
#'   that can be used to estimate the gene-wise dispersion. Common dispersion is
#'   estimated with \code{"optimize"}.
#'   
#'   Arguments that are used by all the methods are:
#'   
#'   \itemize{ \item \code{disp_adjust} \item \code{prop_mode}: Both 
#'   \code{"constrOptim"} and \code{"constrOptimG"} use
#'   \code{\link{constrOptim}} function to maximize the likelihood of
#'   Dirichlet-multinomial proportions. The difference lays in the way the
#'   likelihood and score are computed. \code{"constrOptim"} uses the likelihood
#'   and score that are calculated based on the fact that x*Gamma(x) =
#'   Gamma(x+1). In \code{"constrOptimG"}, we compute them using
#'   \code{\link{lgamma}} function. We recommend using the second approach,
#'   since it is much faster than the first one. \item \code{prop_tol}: The
#'   accuracy for proportions estimation defined as \code{reltol} in
#'   \code{\link{constrOptim}}. }
#'   
#'   Only some of the rest of dispersion parameters in dmDispersion have an 
#'   influence on the output for a given \code{disp_mode}. Here is a list of
#'   such active parameters for different modes:
#'   
#'   \code{"optimize"}, which uses \code{\link{optimize}} to maximize the
#'   profile likelihood.
#'   
#'   \itemize{ \item \code{disp_interval}: Passed as \code{interval}. \item 
#'   \code{disp_tol}: The accuracy defined as \code{tol}. }
#'   
#'   \code{"optim"}, which uses \code{\link{optim}} to maximize the profile 
#'   likelihood.
#'   
#'   \itemize{ \item \code{disp_init} and \code{disp_init_weirMoM}: The initial 
#'   value \code{par}. \item \code{disp_tol}: The accuracy defined as 
#'   \code{factr}. }
#'   
#'   \code{"constrOptim"}, which uses \code{\link{constrOptim}} to maximize the 
#'   profile likelihood.
#'   
#'   \itemize{ \item \code{disp_init} and \code{disp_init_weirMoM}: The initial 
#'   value \code{theta}.. \item \code{disp_tol}: The accuracy defined as 
#'   \code{reltol}. }
#'   
#'   \code{"grid"}, which uses the grid approach from \code{\link{edgeR}}.
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
#'   
#' @param mean_expression Logical. Whether to estimate the mean expression of 
#'   genes.
#' @param common_dispersion Logical. Whether to estimate the common dispersion.
#' @param genewise_dispersion Logical. Whether to estimate the gene-wise 
#'   dispersion.
#' @param disp_adjust Logical. Whether to use the Cox-Reid adjusted or 
#'   non-adjusted profile likelihood.
#' @param disp_mode Optimization method used to maximize the profile likelihood.
#'   Possible values are \code{"optimize"}, \code{"optim"}, 
#'   \code{"constrOptim"}, \code{"grid"}. See Details.
#' @param  disp_interval Numeric vector of length 2 defining the interval of 
#'   possible values for the dispersion.
#' @param disp_tol The desired accuracy when estimating dispersion.
#' @param disp_init Initial dispersion. If \code{common_dispersion} is 
#'   \code{TRUE}, then \code{disp_init} is overwritten by common dispersion 
#'   estimate.
#' @param disp_init_weirMoM Logical. Whether to use the Weir moment estimator as
#'   an initial value for dispersion. If \code{TRUE}, then \code{disp_init} is 
#'   replaced by Weir estimates.
#' @param  disp_grid_length Length of the search grid.
#' @param  disp_grid_range Vector giving the limits of grid interval.
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
#'   values \code{"constrOptim"} and \code{"constrOptimG"}.
#' @param prop_tol The desired accuracy when estimating proportions.
#' @param verbose Numeric. Definie the level of progress messages displayed. 0 -
#'   no messages, 1 - main messages, 2 - message for every gene fitting.
#' @param BPPARAM Parallelization method used by 
#'   \code{\link[BiocParallel]{bplapply}}.
#'   
#' @return Returns a \code{\linkS4class{dmDSdispersion}} or 
#'   \code{\linkS4class{dmSQTLdispersion}} object.
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
#' plotDispersion(d)
#' 
#' head(mean_expression(d))
#' common_dispersion(d)
#' head(genewise_dispersion(d))
#' 
#' }
#' 
#' @seealso \code{\link{data_dmDSdata}}, \code{\link{data_dmSQTLdata}}, 
#'   \code{\link{plotDispersion}}, \code{\link{dmFit}}, \code{\link{dmTest}}
#' @author Malgorzata Nowicka
#' @rdname dmDispersion
#' @export
setMethod("dmDispersion", "dmDSdata", function(x, mean_expression = TRUE, 
  common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = TRUE, 
  disp_mode = "grid", disp_interval = c(0, 1e+5), disp_tol = 1e-08, 
  disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, 
  disp_grid_range = c(-10, 10), disp_moderation = "trended", disp_prior_df = 0, 
  disp_span = 0.1, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = 0, 
  BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  ### Parameter checks:
  stopifnot(is.logical(mean_expression))
  stopifnot(is.logical(common_dispersion))
  stopifnot(is.logical(genewise_dispersion))
  stopifnot(is.logical(disp_adjust))
  stopifnot(length(disp_mode) == 1)
  stopifnot(disp_mode %in% c("optimize", "optim", "constrOptim", "grid"))
  stopifnot(length(disp_interval) == 2)
  stopifnot(disp_interval[1] < disp_interval[2])
  stopifnot(length(disp_tol) == 1)
  stopifnot(is.numeric(disp_tol) && disp_tol > 0)
  stopifnot(length(disp_init) == 1)
  stopifnot(is.numeric(disp_init))
  stopifnot(is.logical(disp_init_weirMoM))
  stopifnot(disp_grid_length > 2)
  stopifnot(length(disp_grid_range) == 2)
  stopifnot(disp_grid_range[1] < disp_grid_range[2])
  stopifnot(length(disp_moderation) == 1)
  stopifnot(disp_moderation %in% c("none", "common", "trended"))
  stopifnot(length(disp_prior_df) == 1)
  stopifnot(is.numeric(disp_prior_df) && disp_prior_df >= 0)
  stopifnot(length(disp_span) == 1)
  stopifnot(is.numeric(disp_span) && disp_span > 0 && disp_span < 1)
  stopifnot(length(prop_mode) == 1)
  stopifnot(prop_mode %in% c("constrOptimG", "constrOptim"))
  stopifnot(length(prop_tol) == 1)
  stopifnot(is.numeric(prop_tol) && prop_tol > 0)
  stopifnot(verbose %in% 0:2)
  
  
  if(mean_expression || (genewise_dispersion && disp_mode == "grid" && 
      disp_moderation == "trended")){
    mean_expression <- dm_estimateMeanExpression(counts = x@counts, 
      verbose = verbose)
  }else{
    mean_expression <- numeric()
  }
  
  if(common_dispersion){
    common_dispersion <- dmDS_estimateCommonDispersion(counts = x@counts, 
      samples = x@samples, disp_adjust = disp_adjust, 
      disp_interval = disp_interval, 
      disp_tol = 1e+01, prop_mode = prop_mode, prop_tol = prop_tol, 
      verbose = verbose, BPPARAM = BPPARAM)
  }else{
    common_dispersion <- numeric()
  }
  
  
  if(genewise_dispersion){
    
    if(length(common_dispersion)){
      message("! Using common_dispersion = ", round(common_dispersion, 2), 
        " as disp_init !")
      disp_init <- common_dispersion
    }
    
    genewise_dispersion <- dmDS_estimateTagwiseDispersion(counts = x@counts, 
      samples = x@samples, mean_expression = mean_expression, 
      disp_adjust = disp_adjust, disp_mode = disp_mode, 
      disp_interval = disp_interval, disp_tol = disp_tol, 
      disp_init = disp_init, 
      disp_init_weirMoM = disp_init_weirMoM, 
      disp_grid_length = disp_grid_length, 
      disp_grid_range = disp_grid_range, disp_moderation = disp_moderation, 
      disp_prior_df = disp_prior_df, disp_span = disp_span, 
      prop_mode = prop_mode, 
      prop_tol = prop_tol, verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    genewise_dispersion <- numeric()
  }
  
  
  return(new("dmDSdispersion", mean_expression = mean_expression, 
    common_dispersion = common_dispersion, 
    genewise_dispersion = genewise_dispersion, 
    counts = x@counts, samples = x@samples))
  
  
})


################################################################################
### plotDispersion
################################################################################

#' Dispersion versus mean expression plot
#' 
#' @return Scatterplot of Dirichlet-multinomial gene-wise dispersion versus mean
#' gene expression. Both variables are scaled with log10. One dot in the plot
#' corresponds to a gene.
#' 
#' @param x \code{\linkS4class{dmDSdispersion}} or
#'   \code{\linkS4class{dmSQTLdispersion}} object.
#' @param ... Other parameters that can be defined by methods using this
#'   generic.
#' @export
setGeneric("plotDispersion", function(x, ...) standardGeneric("plotDispersion"))




##############################


#' @inheritParams plotData
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
#' plotDispersion(d)
#' 
#' }
#' 
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmDSdata}}, \code{\link{data_dmSQTLdata}},
#'   \code{\link{plotData}}, \code{\link{plotFit}}, \code{\link{plotTest}}
#'   
#' @rdname plotDispersion
#' @export
#' @importFrom grDevices pdf dev.off
setMethod("plotDispersion", "dmDSdispersion", function(x, out_dir = NULL){
  
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
  
  if(!is.null(out_dir)){
    pdf(paste0(out_dir, "dispersion_vs_mean.pdf"))
    print(ggp)
    dev.off()
  }else{
    return(ggp)  
  }
  
})








