#' @include class_dmSQTLdispersion.R class_dmDSfit.R
NULL

################################################################################
### dmSQTLfit class
################################################################################

#' dmSQTLfit object
#' 
#' dmSQTLfit extends the \code{\linkS4class{dmDSdispersion}} class by adding the
#' full model Dirichlet-multinomial feature proportion estimates needed for the
#' sQTL analysis. Feature ratios are estimated for each gene and each group that
#' is defined by different SNPs/blocks. Result of \code{\link{dmFit}}.
#' 
#' @slot fit_full List of \code{\linkS4class{MatrixList}} objects. Each element
#'   of this list contains the full model proportion estimates for all the
#'   blocks associated with a given gene.
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmSQTLdata}}, \code{\linkS4class{dmSQTLdata}},
#'   \code{\linkS4class{dmSQTLdispersion}}, \code{\linkS4class{dmSQTLtest}}
setClass("dmSQTLfit", 
  contains = "dmSQTLdispersion",
  representation(fit_full = "list",
    lik_full = "list",
    coef_full = "list",
    fit_full_bb = "list",
    lik_full_bb = "list",
    coef_full_bb = "list"))

########################################


setValidity("dmSQTLfit", function(object){
  # Has to return TRUE when valid object
  
  # TODO: Add checks for other slots
  
  if(!length(object@counts) == length(object@lik_full))
    return("Different number of genes in 'counts' and 'lik_full'")
  
  return(TRUE)
  
})


################################################################################
### show methods
################################################################################

setMethod("show", "dmSQTLfit", function(object){
  
  callNextMethod(object)
  
})


################################################################################
### dmFit
################################################################################


#' @rdname dmFit
#' @export
setMethod("dmFit", "dmSQTLdispersion", function(x, one_way = TRUE, 
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
  
  stopifnot(verbose %in% 0:3)
  
  fit <- dmSQTL_fit(counts = x@counts, genotypes = x@genotypes, 
    dispersion = x@genewise_dispersion,
    one_way = one_way, group_formula = ~ group,
    prop_mode = prop_mode, prop_tol = prop_tol, 
    coef_mode = coef_mode, coef_tol = coef_tol,
    return_fit = FALSE, return_coef = FALSE,
    verbose = verbose, BPPARAM = BPPARAM)
  
  return(new("dmSQTLfit", lik_full = fit[["lik"]], fit_full = fit[["fit"]],
    mean_expression = x@mean_expression, 
    common_dispersion = x@common_dispersion, 
    genewise_dispersion = x@genewise_dispersion, 
    counts = x@counts, genotypes = x@genotypes,
    blocks = x@blocks, samples = x@samples))
  
})


################################################################################
### plotProportions
################################################################################


#' @param snp_id Character indicating the ID of a SNP to be plotted.
#' @rdname plotProportions
#' @export
setMethod("plotProportions", "dmSQTLfit", function(x, gene_id, snp_id, 
  plot_type = "boxplot1", order = TRUE, plot_fit = TRUE, 
  plot_main = TRUE, group_colors = NULL, feature_colors = NULL){
  
  stopifnot(gene_id %in% names(x@blocks))
  
  if(!snp_id %in% x@blocks[[gene_id, "snp_id"]])
    stop(paste0("gene ",gene_id, " and SNP ", snp_id, " do not match!"))
  
  stopifnot(plot_type %in% c("barplot", "boxplot1", "boxplot2", "lineplot", 
    "ribbonplot"))
  stopifnot(is.logical(order))
  stopifnot(is.logical(plot_fit))
  stopifnot(is.logical(plot_main))
  
  counts_gene <- x@counts[[gene_id]]
  block_id <- x@blocks[[gene_id]][x@blocks[[gene_id]][, "snp_id"] == snp_id,
    "block_id"]
  group <- x@genotypes[[gene_id]][block_id, ] 
  
  if(!is.null(group_colors) && 
      plot_type %in% c("barplot", "boxplot1", "lineplot"))
    stopifnot(length(group_colors) == nlevels(group))
  if(!is.null(feature_colors) && 
      plot_type %in% c("boxplot2", "ribbonplot"))
    stopifnot(length(feature_colors) == nrow(counts_gene))
  
  if(nrow(counts_gene) <= 1)
    stop("!Gene has to have at least 2 features! \n")
  
  # Remove NAs
  nonNAs <- !(is.na(counts_gene[1,]) | is.na(group))
  counts_gene <- counts_gene[, nonNAs, drop = FALSE]
  group <- factor(group[nonNAs])
  
  # Order samples by group
  o <- order(group) 
  group <- group[o]
  counts_gene <- counts_gene[, o, drop = FALSE]
  
  main <- NULL
  
  if(plot_main){
    
    mean_expression_gene <- mean(colSums(counts_gene), na.rm = TRUE)
    
    main <- paste0(gene_id, " : ", snp_id, " : ", block_id,
      "\n Mean expression = ", round(mean_expression_gene))
    
    dispersion_gene <- x@genewise_dispersion[[gene_id]][block_id]
    
    main <- paste0(main, ", Precision = ", round(dispersion_gene, 2))
    
  }
  
  prop_full <- NULL
  
  if(plot_fit && length(x@fit_full) > 0){

    fit_full <- x@fit_full[[gene_id]][[which(rownames(x@genotypes[[gene_id]]) == 
        block_id)]][, nonNAs, drop = FALSE]

    prop_full <- fit_full[, !duplicated(group), drop = FALSE]
    colnames(prop_full) <- levels(group)

  }
  
  ggp <- dm_plotProportions(counts = counts_gene, group = group, 
    prop_full = prop_full, main = main, plot_type = plot_type, 
    order = order, group_colors = group_colors, feature_colors = feature_colors)
  
  return(ggp)  
  
})





