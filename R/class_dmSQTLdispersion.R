#' @include class_dmSQTLdata.R
NULL

###############################################################################
### dmSQTLdispersion class
###############################################################################

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
### accessing methods
################################################################################

#' @rdname dmSQTLdispersion-class
#' @export
setMethod("mean_expression", "dmSQTLdispersion", function(x){
  
  data.frame(gene_id = names(x@mean_expression), 
    mean_expression = x@mean_expression, 
    stringsAsFactors = FALSE, row.names = NULL)
  
})


#' @rdname dmSQTLdispersion-class
#' @export
setMethod("common_dispersion", "dmSQTLdispersion", function(x) 
  x@common_dispersion )


#' @rdname dmSQTLdispersion-class
#' @export
setMethod("genewise_dispersion", "dmSQTLdispersion", function(x){
  
  data.frame(gene_id = rep(names(x@genewise_dispersion), 
    sapply(x@genewise_dispersion, length)), 
    block_id = unlist(lapply(x@genewise_dispersion, names)),
    genewise_dispersion = unlist(x@genewise_dispersion), 
    stringsAsFactors = FALSE, row.names = NULL)
  
})

################################################################################
### show methods
################################################################################

setMethod("show", "dmSQTLdispersion", function(object){
  
  callNextMethod(object)
  
  cat("  mean_expression(), common_dispersion(), genewise_dispersion()\n")
  
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
  common_dispersion = TRUE, genewise_dispersion = TRUE, 
  disp_adjust = TRUE, disp_subset = 0.1,
  disp_interval = c(0, 1e+5), disp_tol = 1e+01, 
  disp_init = 100, disp_grid_length = 21, disp_grid_range = c(-10, 10),
  disp_moderation = "none", disp_prior_df = 0, disp_span = 0.1, 
  one_way = TRUE, speed = TRUE,
  prop_mode = "constrOptim", prop_tol = 1e-12, 
  coef_mode = "optim", coef_tol = 1e-12,
  verbose = 0, BPPARAM = BiocParallel::SerialParam()){
  
  ### Parameter checks:
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
    
    ### Use only one SNP per gene and a null model to make computation faster
    genotypes_null <- new("MatrixList", 
      unlistData = matrix(1, nrow = length(genes2keep), 
        ncol = ncol(x@genotypes)), 
      partitioning = split(1:length(genes2keep), 
        factor(names(x@genotypes[genes2keep, ]), 
          levels = names(x@genotypes[genes2keep, ]))))
    
    common_dispersion <- dmSQTL_estimateCommonDispersion(
      counts = x@counts[genes2keep, ], genotypes = genotypes_null, 
      disp_adjust = disp_adjust, 
      disp_interval = disp_interval, disp_tol = disp_tol,
      one_way = one_way, group_formula = ~ 1,
      prop_mode = prop_mode, prop_tol = prop_tol, 
      coef_mode = coef_mode, coef_tol = coef_tol,
      verbose = verbose, BPPARAM = BPPARAM)
    
  }else{
    common_dispersion <- numeric()
  }
  
  
  if(genewise_dispersion){
    
    if(length(common_dispersion)){
      message("! Using common_dispersion = ", round(common_dispersion, 4), 
        " as disp_init !")
      disp_init <- common_dispersion
    }
    
    if(speed){
      
      ### Use only one SNP per gene and a null model to make computation faster
      G <- length(x@genotypes)
      inds <- 1:G
      
      genotypes_null <- new( "MatrixList", 
        unlistData = matrix(1, nrow = G, ncol = ncol(x@genotypes)), 
        partitioning = split(inds, factor(names(x@genotypes), 
          levels = names(x@genotypes))) )
      
      genewise_dispersion <- dmSQTL_estimateTagwiseDispersion(counts = x@counts, 
        genotypes = genotypes_null, mean_expression = mean_expression, 
        disp_adjust = disp_adjust, disp_init = disp_init, 
        disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, 
        disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, 
        disp_span = disp_span, one_way = one_way, group_formula = ~ 1,
        prop_mode = prop_mode, prop_tol = prop_tol, 
        coef_mode = coef_mode, coef_tol = coef_tol,
        verbose = verbose, BPPARAM = BPPARAM)
      
      ### Replicate the values for all the snps
      genewise_dispersion <- relist(rep(unlist(genewise_dispersion), 
        times = elementNROWS(x@genotypes)), x@genotypes@partitioning)
      
    }else{
      
      genewise_dispersion <- dmSQTL_estimateTagwiseDispersion(counts = x@counts, 
        genotypes = x@genotypes, mean_expression = mean_expression, 
        disp_adjust = disp_adjust, disp_init = disp_init, 
        disp_grid_length = disp_grid_length, disp_grid_range = disp_grid_range, 
        disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, 
        disp_span = disp_span, one_way = one_way, group_formula = ~ group,
        prop_mode = prop_mode, prop_tol = prop_tol, 
        coef_mode = coef_mode, coef_tol = coef_tol,
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


###############################################################################
### plotDispersion
###############################################################################


#' @rdname plotDispersion
#' @export
setMethod("plotDispersion", "dmSQTLdispersion", function(x){
  
  if(!length(x@genewise_dispersion) == length(x@counts))
    stop("Genewise dispersion must be estimated for each gene!")
  if(!length(x@genewise_dispersion) == length(x@mean_expression))
    stop("Mean expression must be estimated for each gene!")
  
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
  
  return(ggp)
  
})










































