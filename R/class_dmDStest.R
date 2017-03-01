#' @include class_dmDSfit.R
NULL

###############################################################################
### dmDStest class
###############################################################################

#' dmDStest object
#' 
#' dmDStest extends the \code{\linkS4class{dmDSfit}} class by adding the null 
#' model Dirichlet-multinomial (DM) and beta-binomial (BB) likelihoods and the 
#' gene-level and feature-level results of testing for differential 
#' exon/transcript usage. Result of calling the \code{\link{dmTest}} function.
#' 
#' @return
#' 
#' \itemize{ \item \code{results(x)}: get a data frame with gene-level or 
#' feature-level results.}
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
#' # --------------------------------------------------------------------------
#' # Create dmDSdata object 
#' # --------------------------------------------------------------------------
#' ## Get kallisto transcript counts from the 'PasillaTranscriptExpr' package
#' 
#' library(PasillaTranscriptExpr)
#' 
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
#' 
#' ## Fit full model proportions
#' d <- dmFit(d, design = design)
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
#' 
#' @author Malgorzata Nowicka
#' @seealso \code{\linkS4class{dmDSdata}}, \code{\linkS4class{dmDSdispersion}},
#'   \code{\linkS4class{dmDSfit}}
setClass("dmDStest", 
  contains = "dmDSfit",
  representation(design_fit_null = "matrix",
    lik_null = "numeric",
    lik_null_bb = "numeric",
    results_gene = "data.frame",
    results_feature = "data.frame"))


# ------------------------------------------------------------------------------

setValidity("dmDStest", function(object){
  # Has to return TRUE when valid object!
  
  if(!length(object@lik_null) == length(object@counts))
    return("Different number of genes in 'counts' and 'lik_null'")
  
  if(length(object@lik_null_bb) > 0){
    if(!length(object@lik_null_bb) == nrow(object@counts))
      return("Different number of features in 'counts' and 'lik_null_bb'")
  }
  
  # TODO: Add more checks for results
  
  return(TRUE)
  
})

###############################################################################
### accessing methods
###############################################################################


#' @rdname dmDStest-class
#' @export
setGeneric("results", function(x, ...) standardGeneric("results"))

#' @rdname dmDStest-class
#' @inheritParams dmDSfit-class
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

#' Likelihood ratio test to detect differential transcript/exon usage
#' 
#' First, estimate the null Dirichlet-multinomial and beta-binomial model 
#' parameters and likelihoods using the null model design. Second, perform the 
#' gene-level (DM model) and feature-level (BB model) likelihood ratio tests. In
#' the differential exon/transcript usage analysis, the null model is defined by
#' the null design matrix. In the exon/transcript usage QTL analysis, null
#' models are defined by a design with intercept only. Currently, beta-binomial
#' model is implemented only in the differential usage analysis.
#' 
#' @param x \code{\linkS4class{dmDSfit}} or \code{\linkS4class{dmSQTLfit}} 
#'   object.
#' @param ... Other parameters that can be defined by methods using this 
#'   generic.
#' @export
setGeneric("dmTest", function(x, ...) standardGeneric("dmTest"))


# -----------------------------------------------------------------------------


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
#'   When \code{contrast} is used to define the null model, the null design
#'   matrix is recalculated using the same approach as in
#'   \code{\link[edgeR]{glmLRT}} function from \code{\link{edgeR}}.
#'   
#' @return Returns a \code{\linkS4class{dmDStest}} or 
#'   \code{\linkS4class{dmSQTLtest}} object.
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
#' 
#' ## Fit full model proportions
#' d <- dmFit(d, design = design)
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
#' @seealso \code{\link{plotPValues}} \code{\link[edgeR]{glmLRT}}
#' @references McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression
#' analysis of multifactor RNA-Seq experiments with respect to biological
#' variation. Nucleic Acids Research 40, 4288-4297.
#' @rdname dmTest
#' @export
setMethod("dmTest", "dmDSfit", function(x, 
  coef = NULL, design = NULL, contrast = NULL, 
  one_way = TRUE, bb_model = TRUE,
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
      if(!all(coef %in% colnames(x@design_fit_full))) 
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
  if(bb_model && length(x@lik_full_bb) > 0){
    
    fit0_bb <- bbDS_fit(counts = x@counts, fit = fit0[["fit"]], 
      design = design0, dispersion = x@genewise_dispersion,
      one_way = one_way, verbose = verbose, BPPARAM = BPPARAM)
    
    # Calculate the BB degrees of freedom for the LR test
    df <- rep.int(ncol(x@design_fit_full) - ncol(design0), 
      length(x@lik_full_bb))
    
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
    
  }else{
    
    if(bb_model && length(x@lik_full_bb) == 0)
      message("Beta-Binomial model is not fitted 
        because bb_model=FALSE in dmFit! Rerun dmFit with bb_model=TRUE.")
    
    return(new("dmDStest", 
      results_gene = results_gene,
      design_fit_null = design0, 
      lik_null = fit0[["lik"]],
      design_fit_full = x@design_fit_full, 
      fit_full = x@fit_full, lik_full = x@lik_full, coef_full = x@coef_full,
      lik_full_bb = x@lik_full_bb,  coef_full_bb = x@coef_full_bb,
      mean_expression = x@mean_expression, 
      common_dispersion = x@common_dispersion, 
      genewise_dispersion = x@genewise_dispersion, 
      design_dispersion = x@design_dispersion,
      counts = x@counts, samples = x@samples))
    
  }
  
  
})


###############################################################################
### plotPValues
###############################################################################

#' Plot p-value distribution
#' 
#' @return Plot a histogram of p-values.
#' 
#' @param x \code{\linkS4class{dmDStest}} or \code{\linkS4class{dmSQTLtest}}
#'   object.
#' @export
setGeneric("plotPValues", function(x, ...) standardGeneric("plotPValues"))



# ----------------------------------------------------------------------------




#' @inheritParams results
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
#' 
#' ## Fit full model proportions
#' d <- dmFit(d, design = design)
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
#' @seealso \code{\link{plotData}}, \code{\link{plotDispersion}},
#'   \code{\link{plotProportions}}
#' @rdname plotPValues
#' @export
setMethod("plotPValues", "dmDStest", function(x, level = "gene"){
  
  stopifnot(length(level) == 1)
  stopifnot(level %in% c("gene", "feature"))
  
  res <- slot(x, paste0("results_", level))
  
  if(nrow(res) > 0)
    ggp <- dm_plotPValues(pvalues = res[, "pvalue"])
  else
    stop("Feature-level results are not available! Set bb_model=TRUE in 
      dmFit and dmTest")
  
  return(ggp)  
  
})








































