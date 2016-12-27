#' @include class_dmSQTLdata.R
NULL

################################################################################
### dmSQTLdispersion class
################################################################################

#' dmSQTLdispersion object
#' 
#' dmSQTLdispersion extends the \code{\linkS4class{dmSQTLdata}} by adding the
#' dispersion estimates of Dirichlet-multinomial distribution used to model the
#' feature (e.g., transcript, exon, exonic bin) counts for each gene-SNP pair in
#' the sQTL analysis. Result of \code{\link{dmDispersion}}.
#' 
#' @slot mean_expression Numeric vector of mean gene expression.
#' @slot common_dispersion Numeric value of estimated common dispersion.
#' @slot genewise_dispersion List of estimated gene-wise dispersions. Each
#'   element of this list is a vector of dispersions estimated for all the
#'   genotype blocks assigned to a given gene.
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmSQTLdata}}, \code{\linkS4class{dmSQTLdata}},
#'   \code{\linkS4class{dmSQTLfit}}, \code{\linkS4class{dmSQTLtest}}
setClass("dmSQTLdispersion", 
  contains = "dmSQTLdata",
  representation(mean_expression = "numeric", 
    common_dispersion = "numeric",
    genewise_dispersion = "list"))


#################################


setValidity("dmSQTLdispersion", function(object){
  # has to return TRUE when valid object!
  
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
  
  if(length(object@genewise_dispersion) > 0){
    if(length(object@genewise_dispersion) == length(object@counts)){
      if(all(lapply(object@genewise_dispersion, length) == 
          elementNROWS(object@genotypes)))
        out <- TRUE
      else
        return("Different numbers of blocks in 'genotypes' and in 
          'genewise_dispersion'")
    }
    else 
      return("Unequal number of genes in 'counts' and in 'genewise_dispersion'")
  }
  
  if(length(object@common_dispersion) > 0){
    if(length(object@common_dispersion) == 1)
      out <- TRUE
    else
      return("'common_dispersion' must be a vector of length 1")
  }
  
  return(out)
  
})


################################################################################
### show methods
################################################################################


setMethod("show", "dmSQTLdispersion", function(object){
  
  callNextMethod(object)
  
})


################################################################################
### dmDispersion
################################################################################

#' @rdname dmDispersion
#' @param speed Logical. If \code{FALSE}, dispersion is calculated per each
#'   gene-block. Such calculation may take a long time, since there can be
#'   hundreds of SNPs/blocks per gene. If \code{TRUE}, there will be only one
#'   dipsersion calculated per gene and it will be assigned to all the blocks
#'   matched with this gene.
#' @export
setMethod("dmDispersion", "dmSQTLdata", function(x, mean_expression = TRUE, 
  common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = TRUE, 
  disp_mode = "grid", disp_interval = c(0, 1e+4), disp_tol = 1e-08, 
  disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, 
  disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 0, 
  disp_span = 0.1, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = 0, 
  speed = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  
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
  stopifnot(is.logical(speed))
  
  if(mean_expression || (genewise_dispersion && disp_mode == "grid" && 
      disp_moderation == "trended")){
    mean_expression <- dm_estimateMeanExpression(counts = x@counts, 
      verbose = verbose)
  }else{
    mean_expression <- numeric()
  }
  
  if(common_dispersion){
    
    ### only one SNP per gene (null model)
    inds <- 1:length(x@genotypes)
    genotypes <- new( "MatrixList", 
      unlistData = matrix(1, nrow = length(x@genotypes), ncol = ncol(x@genotypes)), 
      partitioning = split(inds, factor(names(x@genotypes), 
        levels = names(x@genotypes))) )
    
    common_dispersion <- dmSQTL_estimateCommonDispersion(counts = x@counts, 
      genotypes = genotypes, disp_adjust = disp_adjust, disp_interval = disp_interval, 
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
    
    if(speed){
      
      ### only one SNP per gene (null model)
      inds <- 1:length(x@genotypes)
      genotypes <- new( "MatrixList", 
        unlistData = matrix(1, nrow = length(x@genotypes), ncol = ncol(x@genotypes)), 
        partitioning = split(inds, factor(names(x@genotypes), 
          levels = names(x@genotypes))) )
      
      genewise_dispersion <- dmSQTL_estimateTagwiseDispersion(counts = x@counts, 
        genotypes = genotypes, mean_expression = mean_expression, 
        disp_adjust = disp_adjust, disp_mode = disp_mode, 
        disp_interval = disp_interval, disp_tol = disp_tol, 
        disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, 
        disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, 
        disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, 
        disp_span = disp_span, prop_mode = prop_mode, prop_tol = prop_tol, 
        verbose = verbose, BPPARAM = BPPARAM)
      
      ### because we keep only one SNP per gene (null model)
      genewise_dispersion <- relist(rep(unlist(genewise_dispersion), 
        times = elementNROWS(x@genotypes)), x@genotypes@partitioning)
      
    }else{
      
      genewise_dispersion <- dmSQTL_estimateTagwiseDispersion(counts = x@counts, 
        genotypes = x@genotypes, mean_expression = mean_expression, 
        disp_adjust = disp_adjust, disp_mode = disp_mode, 
        disp_interval = disp_interval, disp_tol = disp_tol, 
        disp_init = disp_init, disp_init_weirMoM = disp_init_weirMoM, 
        disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, 
        disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, 
        disp_span = disp_span, prop_mode = prop_mode, prop_tol = prop_tol, 
        verbose = verbose, BPPARAM = BPPARAM)
      
    }
    
  }else{
    genewise_dispersion <- list()
  }
  
  
  return(new("dmSQTLdispersion", mean_expression = mean_expression, 
    common_dispersion = common_dispersion, 
    genewise_dispersion = genewise_dispersion, 
    counts = x@counts, genotypes = x@genotypes, blocks = x@blocks, 
    samples = x@samples))
  
  
})


################################################################################
### plotDispersion
################################################################################


#' @rdname plotDispersion
#' @export
#' @importFrom grDevices pdf dev.off
setMethod("plotDispersion", "dmSQTLdispersion", function(x, out_dir = NULL){
  
  w <- sapply(x@genewise_dispersion, length)
  
  mean_expression <- rep(x@mean_expression, w)
  nr_features <- rep(elementNROWS(x@counts), w)
  
  genewise_dispersion <- unlist(x@genewise_dispersion)
  
  if(length(x@common_dispersion) == 0){
    common_dispersion <- NULL
  }else{
    common_dispersion <- x@common_dispersion
  }
  
  
  ggp <- dm_plotDispersion(genewise_dispersion = genewise_dispersion, 
    mean_expression = mean_expression, nr_features = nr_features, 
    common_dispersion = common_dispersion)
  
  if(!is.null(out_dir)){
    pdf(paste0(out_dir, "dispersion_vs_mean.pdf"))
    print(ggp)
    dev.off()
  }else{
    return(ggp)
  }
  
  
})










































