#' @include class_dmSQTLdata.R
NULL

###############################################################################
### dmSQTLprecision class
###############################################################################

#' dmSQTLprecision object
#' 
#' dmSQTLprecision extends the \code{\linkS4class{dmSQTLdata}} by adding the 
#' precision estimates of Dirichlet-multinomial distribution used to model the 
#' feature (e.g., transcript, exon, exonic bin) counts for each gene-SNP pair in
#' the QTL analysis. Result of \code{\link{dmPrecision}}.
#' 
#' @slot mean_expression Numeric vector of mean gene expression.
#' @slot common_precision Numeric value of estimated common precision.
#' @slot genewise_precision List of estimated gene-wise precisions. Each 
#'   element of this list is a vector of precisions estimated for all the 
#'   genotype blocks assigned to a given gene.
#'   
#' @examples 
#' # --------------------------------------------------------------------------
#' # Create dmSQTLdata object
#' # --------------------------------------------------------------------------
#' # Use subsets of data defined in the GeuvadisTranscriptExpr package
#' 
#' library(GeuvadisTranscriptExpr)
#' \donttest{
#' counts <- GeuvadisTranscriptExpr::counts
#' genotypes <- GeuvadisTranscriptExpr::genotypes
#' gene_ranges <- GeuvadisTranscriptExpr::gene_ranges
#' snp_ranges <- GeuvadisTranscriptExpr::snp_ranges
#' 
#' colnames(counts)[c(1,2)] <- c("feature_id", "gene_id")
#' colnames(genotypes)[4] <- "snp_id"
#' samples <- data.frame(sample_id = colnames(counts)[-c(1,2)])
#' 
#' d <- dmSQTLdata(counts = counts, gene_ranges = gene_ranges,  
#'   genotypes = genotypes, snp_ranges = snp_ranges, samples = samples, 
#'   window = 5e3)
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
#' }
#' @author Malgorzata Nowicka
#' @seealso \code{\linkS4class{dmSQTLdata}}, \code{\linkS4class{dmSQTLfit}},
#'   \code{\linkS4class{dmSQTLtest}}
setClass("dmSQTLprecision", 
  contains = "dmSQTLdata",
  representation(mean_expression = "numeric", 
    common_precision = "numeric",
    genewise_precision = "list"))


# -----------------------------------------------------------------------------


setValidity("dmSQTLprecision", function(object){
  # Has to return TRUE when valid object!
  
  out <- TRUE
  
  if(length(object@mean_expression) > 0){
    if(length(object@mean_expression) == length(object@counts)){
      if(all(names(object@mean_expression) == names(object@counts)))
        out <- TRUE
      else
        return("Different names of 'counts' and 'mean_expression'")
    }
    else 
      return("Unequal length of 'counts' and 'mean_expression'")
  }
  
  if(length(object@genewise_precision) > 0){
    if(length(object@genewise_precision) == length(object@counts)){
      if(all(lapply(object@genewise_precision, length) == 
          elementNROWS(object@genotypes)))
        out <- TRUE
      else
        return("Different numbers of blocks in 'genotypes' and in 
          'genewise_precision'")
    }
    else 
      return("Unequal number of genes in 'counts' and in 'genewise_precision'")
  }
  
  if(length(object@common_precision) > 0){
    if(length(object@common_precision) == 1)
      out <- TRUE
    else
      return("'common_precision' must be a vector of length 1")
  }
  
  return(out)
  
})


################################################################################
### accessing methods
################################################################################

#' @rdname dmSQTLprecision-class
#' @export
setMethod("mean_expression", "dmSQTLprecision", function(x){
  
  data.frame(gene_id = names(x@mean_expression), 
    mean_expression = x@mean_expression, 
    stringsAsFactors = FALSE, row.names = NULL)
  
})


#' @rdname dmSQTLprecision-class
#' @export
setMethod("common_precision", "dmSQTLprecision", function(x) 
  x@common_precision )


#' @rdname dmSQTLprecision-class
#' @export
setMethod("genewise_precision", "dmSQTLprecision", function(x){
  
  data.frame(gene_id = rep(names(x@genewise_precision), 
    sapply(x@genewise_precision, length)), 
    block_id = unlist(lapply(x@genewise_precision, names)),
    genewise_precision = unlist(x@genewise_precision), 
    stringsAsFactors = FALSE, row.names = NULL)
  
})

################################################################################
### show methods
################################################################################

setMethod("show", "dmSQTLprecision", function(object){
  
  callNextMethod(object)
  
  cat("  mean_expression(), common_precision(), genewise_precision()\n")
  
})


################################################################################
### dmPrecision
################################################################################

#' @rdname dmPrecision
#' @param speed Logical. If \code{FALSE}, precision is calculated per each
#'   gene-block. Such calculation may take a long time, since there can be
#'   hundreds of SNPs/blocks per gene. If \code{TRUE}, there will be only one
#'   dipsersion calculated per gene and it will be assigned to all the blocks
#'   matched with this gene.
#' @export
setMethod("dmPrecision", "dmSQTLdata", function(x, mean_expression = TRUE, 
  common_precision = TRUE, genewise_precision = TRUE, 
  prec_adjust = TRUE, prec_subset = 0.1,
  prec_interval = c(0, 1e+5), prec_tol = 1e+01, 
  prec_init = 100, prec_grid_length = 21, prec_grid_range = c(-10, 10),
  prec_moderation = "none", prec_prior_df = 0, prec_span = 0.1, 
  one_way = TRUE, speed = TRUE,
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = 0, BPPARAM = BiocParallel::SerialParam()){
  
  ### Parameter checks:
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
  stopifnot(is.logical(speed))
  
  stopifnot(length(prop_mode) == 1)
  stopifnot(prop_mode %in% c("constrOptim"))
  stopifnot(length(prop_tol) == 1)
  stopifnot(is.numeric(prop_tol) && prop_tol > 0)
  
  stopifnot(length(coef_mode) == 1)
  stopifnot(coef_mode %in% c("optim", "nlminb", "Rcgmin"))
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
    
    ### Use only one SNP per gene and a null model to make computation faster
    genotypes_null <- new("MatrixList", 
      unlistData = matrix(1, nrow = length(genes2keep), 
        ncol = ncol(x@genotypes)), 
      partitioning = split(1:length(genes2keep), 
        factor(names(x@genotypes[genes2keep, ]), 
          levels = names(x@genotypes[genes2keep, ]))))
    
    common_precision <- dmSQTL_estimateCommonPrecision(
      counts = x@counts[genes2keep, ], genotypes = genotypes_null, 
      prec_adjust = prec_adjust, 
      prec_interval = prec_interval, prec_tol = prec_tol,
      one_way = one_way, group_formula = ~ 1,
      prop_mode = prop_mode, prop_tol = prop_tol, 
      coef_mode = coef_mode, coef_tol = coef_tol,
      verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    common_precision <- numeric()
  }
  
  if(genewise_precision){
    
    if(length(common_precision)){
      message("! Using common_precision = ", round(common_precision, 4), 
        " as prec_init !")
      prec_init <- common_precision
    }
    
    if(speed){
      
      ### Use only one SNP per gene and a null model to make computation faster
      G <- length(x@genotypes)
      inds <- 1:G
      
      genotypes_null <- new( "MatrixList", 
        unlistData = matrix(1, nrow = G, ncol = ncol(x@genotypes)), 
        partitioning = split(inds, factor(names(x@genotypes), 
          levels = names(x@genotypes))) )
      
      genewise_precision <- dmSQTL_estimateTagwisePrecision(counts = x@counts, 
        genotypes = genotypes_null, mean_expression = mean_expression, 
        prec_adjust = prec_adjust, prec_init = prec_init, 
        prec_grid_length = prec_grid_length, prec_grid_range = prec_grid_range, 
        prec_moderation = prec_moderation, prec_prior_df = prec_prior_df, 
        prec_span = prec_span, one_way = one_way, group_formula = ~ 1,
        prop_mode = prop_mode, prop_tol = prop_tol, 
        coef_mode = coef_mode, coef_tol = coef_tol,
        verbose = verbose, BPPARAM = BPPARAM)
      
      ### Replicate the values for all the snps
      genewise_precision <- relist(rep(unlist(genewise_precision), 
        times = elementNROWS(x@genotypes)), x@genotypes@partitioning)
      
    }else{
      
      genewise_precision <- dmSQTL_estimateTagwisePrecision(counts = x@counts, 
        genotypes = x@genotypes, mean_expression = mean_expression, 
        prec_adjust = prec_adjust, prec_init = prec_init, 
        prec_grid_length = prec_grid_length, prec_grid_range = prec_grid_range, 
        prec_moderation = prec_moderation, prec_prior_df = prec_prior_df, 
        prec_span = prec_span, one_way = one_way, group_formula = ~ group,
        prop_mode = prop_mode, prop_tol = prop_tol, 
        coef_mode = coef_mode, coef_tol = coef_tol,
        verbose = verbose, BPPARAM = BPPARAM)
      
    }
    
  }else{
    genewise_precision <- list()
  }
  
  return(new("dmSQTLprecision", mean_expression = mean_expression, 
    common_precision = common_precision, 
    genewise_precision = genewise_precision, 
    counts = x@counts, genotypes = x@genotypes, blocks = x@blocks, 
    samples = x@samples))
  
})


###############################################################################
### plotPrecision
###############################################################################


#' @rdname plotPrecision
#' @export
setMethod("plotPrecision", "dmSQTLprecision", function(x){
  
  if(!length(x@genewise_precision) == length(x@counts))
    stop("Genewise precision must be estimated for each gene!")
  if(!length(x@genewise_precision) == length(x@mean_expression))
    stop("Mean expression must be estimated for each gene!")
  
  w <- sapply(x@genewise_precision, length)
  
  mean_expression <- rep(x@mean_expression, w)
  nr_features <- rep(elementNROWS(x@counts), w)
  
  genewise_precision <- unlist(x@genewise_precision)
  
  if(length(x@common_precision) == 0){
    common_precision <- NULL
  }else{
    common_precision <- x@common_precision
  }
  
  ggp <- dm_plotPrecision(genewise_precision = genewise_precision, 
    mean_expression = mean_expression, nr_features = nr_features, 
    common_precision = common_precision)
  
  return(ggp)
  
})










































