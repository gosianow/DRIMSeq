# Prepare data for examples and vignette 

setwd("/home/gosia/R/multinomial_project/package_devel/DRIMSeq")

library(DRIMSeq)
library(devtools)


########################################################
# Examples
########################################################

#############################
### Create dmSQTLdata object
#############################

# Use subsets of data defined in GeuvadisTranscriptExpr package
library(GeuvadisTranscriptExpr)

counts <- GeuvadisTranscriptExpr::counts
genotypes <- GeuvadisTranscriptExpr::genotypes
gene_ranges <- GeuvadisTranscriptExpr::gene_ranges
snp_ranges <- GeuvadisTranscriptExpr::snp_ranges

# Make sure that samples in counts and genotypes are in the same order
sample_id <- colnames(counts[, -(1:2)])

d <- dmSQTLdataFromRanges(counts = counts[, sample_id], gene_id = counts$Gene_Symbol, feature_id = counts$TargetID, gene_ranges = gene_ranges, genotypes = genotypes[, sample_id], snp_id = genotypes$snpId, snp_ranges = snp_ranges, sample_id = sample_id, window = 5e3, BPPARAM = BiocParallel::SerialParam())

plotData(d)

data_dmSQTLdata <- d
use_data(data_dmSQTLdata, overwrite = TRUE)

#############################
### sQTL analysis
#############################
# If possible, use BPPARAM = BiocParallel::MulticoreParam() with more workers

d <- data_dmSQTLdata

head(names(d))
length(d)
d[1:10, ]
d[1:10, 1:10]

### Filtering
d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_expr = 5,  min_samps_feature_prop = 0, minor_allele_freq = 5, BPPARAM = BiocParallel::SerialParam())
plotData(d)


### Calculate dispersion
d <- dmDispersion(d, BPPARAM = BiocParallel::SerialParam())
plotDispersion(d)


### Fit full model proportions
d <- dmFit(d, BPPARAM = BiocParallel::SerialParam())


### Fit null model proportions and test for sQTLs
d <- dmTest(d, BPPARAM = BiocParallel::SerialParam())
plotTest(d)

head(results(d))

### Plot feature proportions for top sQTL
res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

gene_id <- res$gene_id[1]
snp_id <- res$snp_id[1]

plotFit(d, gene_id, snp_id)
plotFit(d, gene_id, snp_id, plot_type = "boxplot2", order = FALSE)
plotFit(d, gene_id, snp_id, plot_type = "ribbonplot")









