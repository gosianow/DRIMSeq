#' @import BiocGenerics
#' @import methods
#' @import BiocParallel
#' @import edgeR
#' @import GenomicRanges
NULL




################################################################################
### Data documentation
################################################################################

#'Sample data for differential splicing analysis
#'
#'We use a subset of kallisto transcript counts from 
#'\code{PasillaTranscriptExpr} package.
#'
#'@format \code{data_dmDSdata} is a \code{\linkS4class{dmDSdata}} object. See 
#'  Examples.
#'  
#'@source Brooks AN, Yang L, Duff MO, et al. Conservation of an RNA regulatory 
#'  map between Drosophila and mammals. Genome Res. 2011;21(2):193-202
#'  
#'  \code{PasillaTranscriptExpr} package
#'  
#'@return \code{data_dmDSdata}
#'  
#' @examples 
#' 
#' #############################
#' ### Create dmDSdata object 
#' #############################
#' ### Get kallisto transcript counts from 'PasillaTranscriptExpr' package
#' 
#' library(PasillaTranscriptExpr)
#' 
#' data_dir  <- system.file("extdata", package = "PasillaTranscriptExpr")
#' 
#' # Load metadata
#' metadata <- read.table(file.path(data_dir, "metadata.txt"), header = TRUE, 
#' as.is = TRUE)
#'
#' # Load counts
#' counts <- read.table(file.path(data_dir, "counts.txt"), header = TRUE, 
#' as.is = TRUE)
#'
#' # Create a samples data frame
#' samples <- data.frame(sample_id = metadata$SampleName, 
#' group = metadata$condition)
#'
#' # Create a dmDSdata object
#' d <- dmDSdata(counts = counts, samples = samples)
#' 
#' plotData(d)
#' 
#' # Use a subset of genes, which is defined in the following file
#' gene_id_subset <- readLines(file.path(data_dir, "gene_id_subset.txt"))
#' d <- d[names(d) %in% gene_id_subset, ]
#' 
#' plotData(d)
#' 
#' data_dmDSdata <- d
#' 
#' ###################################
#' ### Differential splicing analysis
#' ###################################
#' # If possible, use BPPARAM = BiocParallel::MulticoreParam() with more workers
#' 
#' d <- data_dmDSdata
#' 
#' head(counts(d))
#' samples(d)
#' head(names(d))
#' length(d)
#' d[1:20, ]
#' d[1:20, 1:3]
#' 
#' ### Filtering
#' # Check what is the minimal number of replicates per condition 
#' table(samples(d)$group)
#' 
#' d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3, 
#'  min_samps_feature_prop = 0)
#' plotData(d)
#' 
#' ### Calculate dispersion
#' d <- dmDispersion(d, BPPARAM = BiocParallel::SerialParam())
#' plotDispersion(d)
#' 
#' head(mean_expression(d))
#' common_dispersion(d)
#' head(genewise_dispersion(d))
#' 
#' ### Fit full model proportions
#' d <- dmFit(d, BPPARAM = BiocParallel::SerialParam())
#' 
#' head(proportions(d))
#' head(statistics(d))
#' 
#' ### Fit null model proportions and test for DS
#' d <- dmTest(d, BPPARAM = BiocParallel::SerialParam())
#' plotTest(d)
#' 
#' head(proportions(d))
#' head(statistics(d))
#' head(results(d))
#' 
#' ### Plot feature proportions for top DS gene
#' res <- results(d)
#' res <- res[order(res$pvalue, decreasing = FALSE), ]
#' 
#' gene_id <- res$gene_id[1]
#' 
#' plotFit(d, gene_id = gene_id)
#' plotFit(d, gene_id = gene_id, plot_type = "lineplot")
#' plotFit(d, gene_id = gene_id, plot_type = "ribbonplot")
#' 
"data_dmDSdata"





#' Sample data for sQTL analysis
#' 
#' A subset of data from GEUVADIS project where 462 RNA-Seq samples from
#' lymphoblastoid cell lines were obtained. The genome sequencing data of the
#' same individuals is provided by the 1000 Genomes Project. The samples in this
#' project come from five populations: CEPH (CEU), Finns (FIN), British (GBR),
#' Toscani (TSI) and Yoruba (YRI). Here, we use a subset of CEPH data from
#' chromosome 19 available in \code{GeuvadisTranscriptExpr} package.
#' 
#' @format \code{data_dmSQTLdata} is a \code{\linkS4class{dmSQTLdata}} object.
#' See Examples.
#' 
#' @source Lappalainen T, Sammeth M, Friedlander MR, et al. Transcriptome and
#' genome sequencing uncovers functional variation in humans. Nature.
#' 2013;501(7468):506-11
#' 
#' \code{GeuvadisTranscriptExpr} package
#' 
#' @return \code{data_dmSQTLdata}
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
#' data_dmSQTLdata <- d
#' 
#' #############################
#' ### sQTL analysis
#' #############################
#' # If possible, use BPPARAM = BiocParallel::MulticoreParam() with more workers
#' 
#' d <- data_dmSQTLdata
#' 
#' head(names(d))
#' length(d)
#' d[1:10, ]
#' d[1:10, 1:10]
#' 
#' ### Filtering
#' d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_expr = 5, 
#'    min_samps_feature_prop = 0, minor_allele_freq = 5, 
#'    BPPARAM = BiocParallel::SerialParam())
#' plotData(d)
#' 
#' 
#' ### Calculate dispersion
#' d <- dmDispersion(d, BPPARAM = BiocParallel::SerialParam())
#' plotDispersion(d)
#' 
#' ### Fit full model proportions
#' d <- dmFit(d, BPPARAM = BiocParallel::SerialParam())
#' 
#' ### Fit null model proportions and test for sQTLs
#' d <- dmTest(d, BPPARAM = BiocParallel::SerialParam())
#' plotTest(d)
#' 
#' head(results(d))
#' 
#' ### Plot feature proportions for top sQTL
#' res <- results(d)
#' res <- res[order(res$pvalue, decreasing = FALSE), ]
#' 
#' gene_id <- res$gene_id[1]
#' snp_id <- res$snp_id[1]
#' 
#' plotFit(d, gene_id, snp_id)
#' plotFit(d, gene_id, snp_id, plot_type = "boxplot2", order = FALSE)
#' plotFit(d, gene_id, snp_id, plot_type = "ribbonplot")
#' 
"data_dmSQTLdata"











