#' @include class_MatrixList.R
NULL

###############################################################################
### dmSQTLdata class
###############################################################################

#' dmSQTLdata object
#' 
#' dmSQTLdata contains genomic feature expression (counts), genotypes and sample
#' information needed for the sQTL analysis. It can be created with function
#' \code{\link{dmSQTLdata}}.
#' 
#' @return
#' 
#' \itemize{ \item \code{names(x)}: Get the gene names. \item \code{length(x)}:
#' Get the number of genes. \item \code{x[i, j]}: Get a subset of dmDSdata
#' object that consists of counts, genotypes and blocks corresponding to genes i
#' and samples j. }
#' 
#' @param x dmSQTLdata object.
#' @param i,j Parameters used for subsetting.
#'   
#' @slot counts \code{\linkS4class{MatrixList}} of expression, in counts, of 
#'   genomic features. Rows correspond to genomic features, such as exons or 
#'   transcripts. Columns correspond to samples. MatrixList is partitioned in a 
#'   way that each of the matrices in a list contains counts for a single gene.
#' @slot genotypes MatrixList of unique genotypes. Rows correspond to blocks, 
#'   columns to samples. Each matrix in this list is a collection of unique 
#'   genotypes that are matched with a given gene.
#' @slot blocks MatrixList with two columns \code{block_id} and \code{snp_id}. 
#'   For each gene, it identifies SNPs with identical genotypes across the 
#'   samples and assigns them to blocks.
#' @slot samples Data frame with information about samples. It contains unique 
#'   sample names \code{sample_id}.
#'   
#' @examples 
#' #############################
#' ### sQTL analysis
#' #############################
#' 
#' d <- data_dmSQTLdata
#' 
#' head(names(d))
#' length(d)
#' d[1:10, ]
#' d[1:10, 1:10]
#' 
#' @author Malgorzata Nowicka
#' @seealso \code{\link{data_dmSQTLdata}}, 
#'   \code{\linkS4class{dmSQTLdispersion}}, \code{\linkS4class{dmSQTLfit}}, 
#'   \code{\linkS4class{dmSQTLtest}}
setClass("dmSQTLdata", 
  representation(counts = "MatrixList", 
    genotypes = "MatrixList", 
    blocks = "MatrixList",
    samples = "data.frame"))


###################################


setValidity("dmSQTLdata", function(object){
  ### Has to return TRUE when valid object!
  
  if(!ncol(object@counts) == ncol(object@genotypes))
    return(paste0("Unequal number of samples in 'counts' and 'genotypes' ", 
      ncol(object@counts), " and ", ncol(object@genotypes)))
  
  ### Mystery: This does not pass
  # if(!all(colnames(object@blocks) %in% c("block_id", "snp_id")))
  #   return(paste0("'blocks' must contain 'block_id' and 'snp_id' variables"))
  
  if(!all(names(object@counts) == names(object@genotypes)))
    return("'genotypes' and 'counts' do not contain the same genes")
  
  if(!all(names(object@blocks) == names(object@genotypes)))
    return("'genotypes' and 'blocks' do not contain the same genes or SNPs")
  
  return(TRUE)
  
})


###############################################################################
### accessing and subsetting methods
###############################################################################


setMethod("show", "dmSQTLdata", function(object){
  
  cat("An object of class", class(object), "\n")
  
  cat("with", length(object), "genes and", ncol(object@counts), "samples\n")
  
})



###############################

#' @rdname dmSQTLdata-class
#' @export
setMethod("names", "dmSQTLdata", function(x) names(x@counts) )

#' @rdname dmSQTLdata-class
#' @export
setMethod("length", "dmSQTLdata", function(x) length(x@counts) )


#' @aliases [,dmSQTLdata-method [,dmSQTLdata,ANY-method
#' @rdname dmSQTLdata-class
#' @export
setMethod("[", "dmSQTLdata", function(x, i, j){
  
  if(missing(j)){
    
    counts <- x@counts[i, , drop = FALSE]
    genotypes <- x@genotypes[i, , drop = FALSE]
    blocks <- x@blocks[i, , drop = FALSE]
    samples <- x@samples
    
  }else{
    
    if(missing(i)){
      counts <- x@counts[, j, drop = FALSE]
      genotypes <- x@genotypes[, j, drop = FALSE]
    }else{
      counts <- x@counts[i, j, drop = FALSE]
      genotypes <- x@genotypes[i, j, drop = FALSE]
      blocks <- x@blocks[i, , drop = FALSE]
    }
    
    samples <- x@samples
    rownames(samples) <- samples$sample_id
    samples <- samples[j, , drop = FALSE]
    samples$sample_id <- factor(samples$sample_id)
    rownames(samples) <- NULL
    
  }
  
  return(new("dmSQTLdata", counts = counts, genotypes = genotypes, 
    blocks = blocks, samples = samples))
  
})

###############################################################################
### dmSQTLdata
###############################################################################

blocks_per_gene <- function(g, genotypes){
  # g = 1
  
  genotypes_df <- data.frame(t(genotypes[[g]]))
  matching_snps <- match(genotypes_df, genotypes_df)
  oo <- order(matching_snps, decreasing = FALSE)
  block_id <- paste0("block_", as.numeric(factor(matching_snps)))
  snp_id <- colnames(genotypes_df)
  blocks_tmp <- cbind(block_id, snp_id)
  
  return(blocks_tmp[oo, , drop = FALSE])
  
}

#' Create dmSQTLdata object
#' 
#' Constructor functions for a \code{\linkS4class{dmSQTLdata}} object. 
#' dmSQTLdata assignes to a gene all the SNPs that are located in a given
#' surrounding (\code{window}) of this gene.
#' 
#' It is quite common that sample grouping defined by some of the SNPs is 
#' identical. Compare \code{dim(genotypes)} and \code{dim(unique(genotypes))}. 
#' In our sQTL analysis, we do not repeat tests for the SNPs that define the 
#' same grouping of samples. Each grouping is tested only once. SNPs that define
#' such unique groupings are aggregated into blocks. P-values and adjusted 
#' p-values are estimated at the block level, but the returned results are 
#' extended to a SNP level by repeating the block statistics for each SNP that 
#' belongs to a given block.
#' 
#' @inheritParams dmDSdata
#' @param genotypes Data frame with genotypes. Rows correspond to SNPs. This 
#'   data frame has to contain a \code{snp_id} column with SNP IDs and columns 
#'   with genotypes for each sample. Column names corresponding to sample IDs 
#'   must be the same as in the \code{sample} data frame. The genotype of each 
#'   sample is coded in the following way: 0 for ref/ref, 1 for ref/not ref, 2 
#'   for not ref/not ref, -1 or \code{NA} for missing value.
#' @param gene_ranges \code{\linkS4class{GRanges}} object with gene location. It
#'   must contain gene names when calling names().
#' @param snp_ranges \code{\linkS4class{GRanges}} object with SNP location. It 
#'   must contain SNP names when calling names().
#' @param window Size of a down and up stream window, which is defining the 
#'   surrounding for a gene. Only SNPs that are located within a gene or its 
#'   surrounding are considered in the sQTL analysis.
#' @param samples Data frame with column \code{sample_id} corresponding to 
#'   unique sample IDs
#' @param BPPARAM Parallelization method used by 
#'   \code{\link[BiocParallel]{bplapply}}.
#'   
#' @return Returns a \code{\linkS4class{dmSQTLdata}} object.
#'   
#' @examples 
#'  
#' #############################
#' ### Create dmSQTLdata object
#' #############################
#' 
#' # Use subsets of data defined in GeuvadisTranscriptExpr package
#' library(GeuvadisTranscriptExpr)
#' 
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
#' genotypes = genotypes, snp_ranges = snp_ranges, samples = samples, 
#' window = 5e3, BPPARAM = BiocParallel::SerialParam())
#' 
#' plotData(d)
#' 
#' @seealso \code{\link{data_dmSQTLdata}}, \code{\link{dmFilter}}, 
#'   \code{\link{dmDispersion}}, \code{\link{dmFit}}, \code{\link{dmTest}}
#' @author Malgorzata Nowicka
#' @export
#' @importFrom IRanges width
#' @importFrom S4Vectors queryHits subjectHits
dmSQTLdata <- function(counts, gene_ranges, genotypes, snp_ranges, samples, 
  window = 5e3, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  stopifnot(is.numeric(window))
  stopifnot(window >= 0)
  
  ### Check on samples
  stopifnot(class(samples) == "data.frame")
  stopifnot("sample_id" %in% colnames(samples))
  stopifnot(sum(duplicated(samples$sample_id)) == 0)
  
  ### Check on counts
  stopifnot(class(counts) == "data.frame")
  stopifnot(all(c("gene_id", "feature_id") %in% colnames(counts)))
  stopifnot(all(samples$sample_id %in% colnames(counts)))
  
  ### Check on genotypes
  stopifnot(class(genotypes) == "data.frame")
  stopifnot("snp_id" %in% colnames(genotypes))
  stopifnot(all(samples$sample_id %in% colnames(counts)))
  
  sample_id <- samples$sample_id
  gene_id <- counts$gene_id
  feature_id <- counts$feature_id
  snp_id <- genotypes$snp_id
  
  stopifnot( class( gene_id ) %in% c("character", "factor"))
  stopifnot( class( feature_id ) %in% c("character", "factor"))
  stopifnot( class( sample_id ) %in% c("character", "factor"))
  stopifnot( class( snp_id ) %in% c("character", "factor"))
  
  stopifnot(all(!is.na(gene_id)))
  stopifnot(all(!is.na(feature_id)))
  stopifnot(all(!is.na(sample_id)))
  stopifnot(all(!is.na(snp_id)))
  
  counts <- counts[, as.character(sample_id), drop = FALSE]
  counts <- as.matrix(counts)
  stopifnot(mode(counts) %in% "numeric")
  
  genotypes <- genotypes[, as.character(sample_id), drop = FALSE]
  genotypes <- as.matrix(genotypes)
  stopifnot(mode(genotypes) %in% "numeric")
  stopifnot(all(genotypes %in% c(-1, 0, 1, 2, NA)))
  rownames(genotypes) <- snp_id
  
  stopifnot(class(gene_ranges) == "GRanges")
  stopifnot(class(snp_ranges) == "GRanges")
  stopifnot(!is.null(names(snp_ranges)))
  stopifnot(!is.null(names(gene_ranges)))
  
  ### Keep genes that are in counts and in gene_ranges
  genes_overlap <- intersect(names(gene_ranges), gene_id)
  
  genes2keep <- gene_id %in% genes_overlap
  counts <- counts[genes2keep, , drop = FALSE]
  gene_id <- gene_id[genes2keep]
  feature_id <- feature_id[genes2keep]
  
  gene_ranges <- gene_ranges[genes_overlap, ]
  
  ### Keep SNPs that are in genotypes and in snp_ranges 
  ### and make them in the same order
  snps_overlap <- intersect(names(snp_ranges), snp_id)
  
  genotypes <- genotypes[snps_overlap, , drop = FALSE]
  snp_ranges <- snp_ranges[snps_overlap, ]
  
  gene_ranges <- GenomicRanges::resize(gene_ranges, 
    width(gene_ranges) + 2 * window, fix = "center")
  
  ## Match genes and SNPs
  variantMatch <- GenomicRanges::findOverlaps(gene_ranges, snp_ranges, 
    select = "all")
  
  q <- queryHits(variantMatch)
  s <- subjectHits(variantMatch)
  
  genotypes <- genotypes[s, ]
  snp_id <- snp_id[s]
  gene_id_genotypes <- names(gene_ranges)[q]
  
  
  ### keep genes that are in counts and in genotypes
  genes2keep <- gene_id %in% gene_id_genotypes
  counts <- counts[genes2keep, , drop = FALSE]
  gene_id <- gene_id[genes2keep]
  feature_id <- feature_id[genes2keep]
  
  genes2keep <- gene_id_genotypes %in% gene_id
  genotypes <- genotypes[genes2keep, , drop = FALSE]
  gene_id_genotypes <- gene_id_genotypes[genes2keep]
  snp_id <- snp_id[genes2keep]
  
  ### order genes in counts and in genotypes
  if(class(gene_id) == "character")
    gene_id <- factor(gene_id, levels = unique(gene_id))
  
  order_counts <- order(gene_id)
  counts <- counts[order_counts, , drop = FALSE]
  gene_id <- gene_id[order_counts]
  feature_id <- feature_id[order_counts]
  
  gene_id_genotypes <- factor(gene_id_genotypes, levels = levels(gene_id))
  order_genotypes <- order(gene_id_genotypes)
  genotypes <- genotypes[order_genotypes, , drop = FALSE]
  gene_id_genotypes <- gene_id_genotypes[order_genotypes]
  snp_id <- snp_id[order_genotypes]
  
  colnames(counts) <- sample_id
  rownames(counts) <- feature_id
  
  colnames(genotypes) <- sample_id
  rownames(genotypes) <- snp_id
  
  inds_counts <- 1:length(gene_id)
  names(inds_counts) <- feature_id
  partitioning_counts <- split(inds_counts, gene_id)
  
  inds_genotypes <- 1:length(gene_id_genotypes)
  names(inds_genotypes) <- snp_id
  partitioning_genotypes <- split(inds_genotypes, gene_id_genotypes)
  
  counts <- new( "MatrixList", unlistData = counts, 
    partitioning = partitioning_counts)
  genotypes <- new( "MatrixList", unlistData = genotypes, 
    partitioning = partitioning_genotypes)
  
  ### Keep unique genotypes and create info about blocs
  inds <- 1:length(genotypes)
  
  blocks <- MatrixList(BiocParallel::bplapply(inds, blocks_per_gene, 
    genotypes = genotypes, BPPARAM = BPPARAM))
  
  names(blocks) <- names(genotypes)
  
  genotypes_u <- MatrixList(lapply(inds, function(g){
    # g = 1
    
    genotypes_tmp <- unique(genotypes[[g]])
    rownames(genotypes_tmp) <- paste0("block_", 1:nrow(genotypes_tmp))
    return(genotypes_tmp)
    
  }))
  
  names(genotypes_u) <- names(genotypes)
  samples <- data.frame(sample_id = sample_id)
  data <- new("dmSQTLdata", counts = counts, genotypes = genotypes_u, 
    blocks = blocks, samples = samples)
  
  return(data)
  
  
}


################################################################################
### dmFilter
################################################################################


#' @param minor_allele_freq Minimal number of samples where each of the
#'   genotypes has to be present.
#' @param BPPARAM Parallelization method used by
#'   \code{\link[BiocParallel]{bplapply}}.
#' @details
#' 
#' In sQTL analysis, usually, we deal with data that has many more replicates
#' than data from a standard differential splicing assay. Our example data set
#' consists of 91 samples. Requiring that genes are expressed in all samples may
#' be too stringent, especially since there may be missing values in the data
#' and for some genes you may not observe counts in all 91 samples. Slightly
#' lower threshold ensures that we do not eliminate such genes. For example, if
#' \code{min_samps_gene_expr = 70} and \code{min_gene_expr = 10}, only genes
#' with expression of at least 10 in at least 70 samples are kept. Samples with
#' expression lower than 10 have \code{NA}s assigned and are skipped in the
#' analysis of this gene. \code{minor_allele_freq} indicates the minimal number
#' of samples for the minor allele presence. Usually, it is equal to 5\% of
#' total samples.
#' 
#' @examples 
#' #############################
#' ### sQTL analysis
#' #############################
#' # If possible, use BPPARAM = BiocParallel::MulticoreParam() with more workers
#' 
#' d <- data_dmSQTLdata
#' \donttest{
#' ### Filtering
#' d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_expr = 5, 
#'  min_samps_feature_prop = 0, minor_allele_freq = 5, 
#'  BPPARAM = BiocParallel::SerilaParam())
#' plotData(d)
#' }
#' @rdname dmFilter
#' @export
setMethod("dmFilter", "dmSQTLdata", function(x, min_samps_gene_expr, 
  min_samps_feature_expr, min_samps_feature_prop, minor_allele_freq, 
  min_gene_expr = 10, min_feature_expr = 10, min_feature_prop = 0, 
  max_features = Inf, BPPARAM = BiocParallel::MulticoreParam(workers = 1)){
  
  stopifnot(min_samps_gene_expr >= 0 && 
      min_samps_gene_expr <= ncol(x@counts))
  stopifnot(min_gene_expr >= 0)
  stopifnot(min_samps_feature_expr >= 0 && 
      min_samps_feature_expr <= ncol(x@counts))
  stopifnot(min_feature_expr >= 0)
  stopifnot(min_samps_feature_prop >= 0 && 
      min_samps_feature_prop <= ncol(x@counts))
  stopifnot(min_feature_prop >= 0 && min_feature_prop <= 1)
  stopifnot(max_features >= 2)
  stopifnot(minor_allele_freq >= 1 && 
      minor_allele_freq <= floor(ncol(x@counts)/2))
  
  data_filtered <- dmSQTL_filter(counts = x@counts, genotypes = x@genotypes, 
    blocks = x@blocks, samples = x@samples, 
    min_samps_gene_expr = min_samps_gene_expr, 
    min_gene_expr = min_gene_expr, 
    min_samps_feature_expr = min_samps_feature_expr, 
    min_feature_expr = min_feature_expr, 
    min_samps_feature_prop = min_samps_feature_prop, 
    min_feature_prop = min_feature_prop, max_features = max_features,
    minor_allele_freq = minor_allele_freq, BPPARAM = BPPARAM)
  
  return(data_filtered)
  
  
})


###############################################################################
### plotData
###############################################################################


#' @examples 
#' 
#' #############################
#' ### sQTL analysis
#' #############################
#' 
#' d <- data_dmSQTLdata
#' plotData(d)
#'
#' @rdname plotData
#' @export
#' @importFrom grDevices pdf dev.off
setMethod("plotData", "dmSQTLdata", function(x, out_dir = NULL){
  
  tt <- elementNROWS(x@counts)
  ggp1 <- dm_plotDataFeatures(tt)
  
  
  tt <- elementNROWS(x@blocks)
  ggp2 <- dm_plotDataSnps(tt)
  
  
  tt <- elementNROWS(x@genotypes)
  ggp3 <- dm_plotDataBlocks(tt)
  
  
  if(!is.null(out_dir)){
    pdf(paste0(out_dir, "hist_features.pdf"))
    print(ggp1)
    dev.off()
    
    pdf(paste0(out_dir, "hist_snps.pdf"))
    print(ggp2)
    dev.off()
    
    pdf(paste0(out_dir, "hist_blocks.pdf"))
    print(ggp3)
    dev.off()
  }else{
    return(list(ggp1, ggp2, ggp3))
  }
  
  
})





























