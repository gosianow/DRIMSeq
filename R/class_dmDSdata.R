#' @include class_MatrixList.R
NULL

###############################################################################
### dmDSdata class
###############################################################################

#' dmDSdata object
#' 
#' dmDSdata contains expression, in counts, of genomic features such as exons or
#' transcripts and sample information needed for the differential 
#' exon/transcript usage (DEU or DTU) analysis. It can be created with function 
#' \code{\link{dmDSdata}}.
#' 
#' @return
#' 
#' \itemize{ \item \code{counts(object)}: Get a data frame with counts. \item 
#' \code{samples(x)}: Get a data frame with the sample information. \item 
#' \code{names(x)}: Get the gene names. \item \code{length(x)}: Get the number 
#' of genes. \item \code{x[i, j]}: Get a subset of dmDSdata object that consists
#' of counts for genes i and samples j. }
#' 
#' @param object,x dmDSdata object.
#' @param i,j Parameters used for subsetting.
#' @param ... Other parameters that can be defined by methods using this 
#'   generic.
#'   
#' @slot counts \code{\linkS4class{MatrixList}} of expression, in counts, of 
#'   genomic features. Rows correspond to genomic features, such as exons or 
#'   transcripts. Columns correspond to samples. MatrixList is partitioned in a 
#'   way that each of the matrices in a list contains counts for a single gene.
#' @slot samples Data frame with information about samples. It must contain 
#'   \code{sample_id} variable with unique sample names and other covariates 
#'   that desribe samples and are needed for the differential analysis.
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
#' }
#' @author Malgorzata Nowicka
#' @seealso \code{\linkS4class{dmDSprecision}}, \code{\linkS4class{dmDSfit}},
#'   \code{\linkS4class{dmDStest}}
setClass("dmDSdata", 
  representation(counts = "MatrixList", samples = "data.frame"))


###############################


setValidity("dmDSdata", function(object){
  # has to return TRUE when valid object!
  
  if(!ncol(object@counts) == nrow(object@samples))
    return(paste0("Unequal number of samples in 'counts' and 'samples' ", 
      ncol(object@counts), " and ", nrow(object@samples), "!"))
  
  if(!all(c("sample_id") %in% colnames(object@samples)))
    return("'samples' must contain 'sample_id' variable!")
  
  if(!length(unique(object@samples$sample_id)) == nrow(object@samples))
    return("There must be a unique 'sample_id' for each sample!")
  
  if(!all(colnames(object@counts) == object@samples$sample_id))
    return("Column names of 'counts' must be the same as 'sample_id' 
      in 'samples'!")
  
  return(TRUE)
  
})


###############################################################################
### show, accessing and subsetting methods
###############################################################################

#' @rdname dmDSdata-class
#' @export
setMethod("counts", "dmDSdata", function(object){
  
  data.frame(gene_id = rep(names(object@counts), elementNROWS(object@counts)), 
    feature_id = rownames(object@counts), 
    object@counts@unlistData, 
    stringsAsFactors = FALSE, row.names = NULL)
  
})


#' @rdname dmDSdata-class
#' @export
setGeneric("samples", function(x, ...) standardGeneric("samples"))


#' @rdname dmDSdata-class
#' @export
setMethod("samples", "dmDSdata", function(x) x@samples )



################################

setMethod("show", "dmDSdata", function(object){
  
  cat("An object of class", class(object), "\n")
  
  cat("with", length(object), "genes and", ncol(object@counts), "samples\n")
  
  cat("* data accessors: counts(), samples()\n")
  
})

################################

#' @rdname dmDSdata-class
#' @export
setMethod("names", "dmDSdata", function(x) names(x@counts) )

#' @rdname dmDSdata-class
#' @export
setMethod("length", "dmDSdata", function(x) length(x@counts) )


#' @aliases [,dmDSdata-method [,dmDSdata,ANY-method
#' @rdname dmDSdata-class
#' @export
setMethod("[", "dmDSdata", function(x, i, j){
  
  if(missing(j)){
    
    counts <- x@counts[i, , drop = FALSE]
    samples <- x@samples
    
  }else{
    
    if(missing(i)){
      counts <- x@counts[, j, drop = FALSE]
    }else{
      counts <- x@counts[i, j, drop = FALSE]
    }
    
    samples <- x@samples
    rownames(samples) <- samples$sample_id
    samples <- samples[j, , drop = FALSE]
    
    # Drop unused levels for factors
    for(i in 1:ncol(samples)){
      if(class(samples[, i]) == "factor")
        samples[, i] <- factor(samples[, i])
    }
    
    rownames(samples) <- NULL
    
  }
  
  return(new("dmDSdata", counts = counts, samples = samples))
  
})


###############################################################################
### dmDSdata
###############################################################################

#' Create dmDSdata object
#' 
#' Constructor function for a \code{\linkS4class{dmDSdata}} object.
#' 
#' @param counts Data frame with counts. Rows correspond to features, for
#'   example, transcripts or exons. This data frame has to contain a
#'   \code{gene_id} column with gene IDs, \code{feature_id} column with feature
#'   IDs and columns with counts for each sample. Column names corresponding to
#'   sample IDs must be the same as in the \code{sample} data frame.
#' @param samples Data frame where each row corresponds to one sample. Columns
#'   have to contain unique sample IDs in \code{sample_id} variable and a
#'   grouping variable \code{group}.
#' @return Returns a \linkS4class{dmDSdata} object.
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
#' }
#' @seealso \code{\link{plotData}}
#' @author Malgorzata Nowicka
#' @export
dmDSdata <- function(counts, samples){
  
  ### Check on samples
  stopifnot(class(samples) == "data.frame")
  stopifnot("sample_id" %in% colnames(samples))
  stopifnot(sum(duplicated(samples$sample_id)) == 0)
  
  ### Check on counts
  stopifnot(class(counts) == "data.frame")
  stopifnot(all(c("gene_id", "feature_id") %in% colnames(counts)))
  stopifnot(all(samples$sample_id %in% colnames(counts)))
  stopifnot(sum(duplicated(counts$feature_id)) == 0)
  
  gene_id <- counts$gene_id
  feature_id <- counts$feature_id
  
  stopifnot(class(gene_id) %in% c("character", "factor"))
  stopifnot(class(feature_id) %in% c("character", "factor"))
  stopifnot(class(samples$sample_id) %in% c("character", "factor"))
  
  stopifnot(all(!is.na(gene_id)))
  stopifnot(all(!is.na(feature_id)))
  stopifnot(all(!is.na(samples$sample_id)))
  
  counts <- counts[, as.character(samples$sample_id), drop = FALSE]
  
  counts <- as.matrix(counts)
  stopifnot(mode(counts) %in% "numeric")
  
  if(class(gene_id) == "character")
    gene_id <- factor(gene_id, levels = unique(gene_id))
  else 
    gene_id <- factor(gene_id)
  
  for(i in 1:ncol(samples)){
    
    if(class(samples[, i]) == "character")
      samples[, i] <- factor(samples[, i], levels = unique(samples[, i]))
    else if(class(samples[, i]) == "factor")
      samples[, i] <- factor(samples[, i])
    
  }
  
  # Ordering
  or <- order(gene_id)
  
  counts <- counts[or, , drop = FALSE]
  gene_id <- gene_id[or]
  feature_id <- feature_id[or]
  
  rownames(counts) <- feature_id
  
  inds <- 1:length(gene_id)
  names(inds) <- feature_id
  
  partitioning <- split(inds, gene_id)
  
  data <- new("dmDSdata", counts = new("MatrixList", unlistData = counts, 
    partitioning = partitioning), samples = samples)
  
  return(data)
  
}


###############################################################################
### dmFilter
###############################################################################



#' Filtering
#' 
#' Filtering of genes and features with low expression. Additionally, for the
#' dmSQTLdata object, filtering of genotypes with low frequency.
#' 
#' @param x \code{\linkS4class{dmDSdata}} or \code{\linkS4class{dmSQTLdata}}
#'   object.
#' @param ... Other parameters that can be defined by methods using this
#'   generic.
#' @export
setGeneric("dmFilter", function(x, ...) standardGeneric("dmFilter"))

# -----------------------------------------------------------------------------

#' @details Filtering parameters should be adjusted according to the sample size
#'   of the experiment data and the number of replicates per condition.
#'   
#'   \code{min_samps_gene_expr} defines the minimal number of samples where 
#'   genes are required to be expressed at the minimal level of 
#'   \code{min_gene_expr} in order to be included in the downstream analysis. 
#'   Ideally, we would like that genes were expressed at some minimal level in 
#'   all samples because this would lead to better estimates of feature ratios.
#'   
#'   Similarly, \code{min_samps_feature_expr} and \code{min_samps_feature_prop} 
#'   defines the minimal number of samples where features are required to be 
#'   expressed at the minimal levels of counts \code{min_feature_expr} or 
#'   proportions \code{min_feature_prop}. In differential transcript/exon usage
#'   analysis, we suggest using \code{min_samps_feature_expr} and 
#'   \code{min_samps_feature_prop} equal to the minimal number of replicates in 
#'   any of the conditions. For example, in an assay with 3 versus 5 replicates,
#'   we would set these parameters to 3, which allows a situation where a 
#'   feature is expressed in one condition but may not be expressed at all in 
#'   another one, which is an example of differential transcript/exon usage.
#'   
#'   By default, all the filtering parameters equal zero which means that 
#'   features with zero expression in all samples are removed as well as genes 
#'   with only one non-zero feature.
#'   
#' @param min_samps_gene_expr Minimal number of samples where genes should be 
#'   expressed. See Details.
#' @param min_gene_expr Minimal gene expression.
#' @param min_samps_feature_expr Minimal number of samples where features should
#'   be expressed. See Details.
#' @param min_feature_expr Minimal feature expression.
#' @param min_samps_feature_prop Minimal number of samples where features should
#'   be expressed. See details.
#' @param min_feature_prop Minimal proportion for feature expression. This value
#'   should be between 0 and 1.
#' @param run_gene_twice Whether to re-run the gene-level filter
#'   after the feature-level filters.
#' @return Returns filtered \code{\linkS4class{dmDSdata}} or 
#'   \code{\linkS4class{dmSQTLdata}} object.
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
#' }
#' @seealso \code{\link{plotData}}
#' @author Malgorzata Nowicka
#' @rdname dmFilter
#' @export
setMethod("dmFilter", "dmDSdata", function(x, min_samps_gene_expr = 0, 
  min_samps_feature_expr = 0, min_samps_feature_prop = 0, min_gene_expr = 0, 
  min_feature_expr = 0, min_feature_prop = 0, run_gene_twice = FALSE){
  
  stopifnot(min_samps_gene_expr >= 0 && 
      min_samps_gene_expr <= ncol(x@counts))
  stopifnot(min_gene_expr >= 0)
  stopifnot(min_samps_feature_expr >= 0 && 
      min_samps_feature_expr <= ncol(x@counts))
  stopifnot(min_feature_expr >= 0)
  stopifnot(min_samps_feature_prop >= 0 && 
      min_samps_feature_prop <= ncol(x@counts))
  stopifnot(min_feature_prop >= 0 && min_feature_prop <= 1)
  
  counts_filtered <- dmDS_filter(counts = x@counts, 
    min_samps_gene_expr = min_samps_gene_expr,  
    min_gene_expr = min_gene_expr, 
    min_samps_feature_expr = min_samps_feature_expr, 
    min_feature_expr = min_feature_expr,
    min_samps_feature_prop = min_samps_feature_prop, 
    min_feature_prop = min_feature_prop,
    run_gene_twice = run_gene_twice)
  
  return(new("dmDSdata", counts = counts_filtered, samples = x@samples))
  
})


################################################################################
### plotData
################################################################################


#' Plot data summary
#' 
#' @return Returns a \code{ggplot} object and can be further modified, for
#'   example, using \code{theme()}. Plots a histogram of the number of features
#'   per gene. Additionally, for \code{\linkS4class{dmSQTLdata}} object, plots a
#'   histogram of the number of SNPs per gene and a histogram of the number of
#'   unique SNPs (blocks) per gene.
#'   
#' @param x \code{\linkS4class{dmDSdata}} or \code{\linkS4class{dmSQTLdata}} 
#'   object.
#' @param ... Other parameters that can be defined by methods using this 
#'   generic.
#' @export
setGeneric("plotData", function(x, ...) standardGeneric("plotData"))


# -----------------------------------------------------------------------------


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
#' }
#' @author Malgorzata Nowicka
#' @seealso \code{\link{plotPrecision}}, \code{\link{plotProportions}},
#'   \code{\link{plotPValues}}
#' @rdname plotData
#' @export
setMethod("plotData", "dmDSdata", function(x){
  
  tt <- elementNROWS(x@counts)
  ggp <- dm_plotDataFeatures(tt = tt)
  
  return(ggp)
  
})



