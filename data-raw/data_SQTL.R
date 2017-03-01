# Prepare data for examples and vignette 

setwd("/home/gosia/R/package_devel/DRIMSeq")

library(DRIMSeq)
library(devtools)

load_all("/home/gosia/R/package_devel/DRIMSeq")


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

colnames(counts)[c(1,2)] <- c("feature_id", "gene_id")
colnames(genotypes)[4] <- "snp_id"
samples <- data.frame(sample_id = colnames(counts)[-c(1,2)])

d <- dmSQTLdata(counts = counts, gene_ranges = gene_ranges,  
  genotypes = genotypes, snp_ranges = snp_ranges, samples = samples, 
  window = 5e3, BPPARAM = BiocParallel::SerialParam())

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
d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_expr = 5,  
  min_samps_feature_prop = 0, minor_allele_freq = 5, 
  BPPARAM = BiocParallel::SerialParam())
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









