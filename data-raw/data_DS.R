# Prepare data for examples and vignette 

setwd("/home/gosia/R/package_devel/DRIMSeq")

library(DRIMSeq)
library(devtools)

load_all("/home/gosia/R/package_devel/DRIMSeq")

########################################################
# Examples
########################################################


#############################
### Create dmDSdata object
#############################
### Get kallisto transcript counts from 'PasillaTranscriptExpr' package

library(PasillaTranscriptExpr)

data_dir  <- system.file("extdata", package = "PasillaTranscriptExpr")

# Load metadata
metadata <- read.table(file.path(data_dir, "metadata.txt"), header = TRUE, 
  as.is = TRUE)

# Load counts
counts <- read.table(file.path(data_dir, "counts.txt"), header = TRUE, 
  as.is = TRUE)

# Create a samples data frame
samples <- data.frame(sample_id = metadata$SampleName, 
  group = metadata$condition)

# Create a dmDSdata object
d <- dmDSdata(counts = counts, samples = samples)
d

plotData(d)

# Use a subset of genes, which is defined in the following file
gene_id_subset <- readLines(file.path(data_dir, "gene_id_subset.txt"))
d <- d[names(d) %in% gene_id_subset, ]
d

plotData(d)

data_dmDSdata <- d
use_data(data_dmDSdata, overwrite = TRUE)


###################################
### Differential splicing analysis
###################################
# If possible, use BPPARAM = BiocParallel::MulticoreParam() with more workers

d <- data_dmDSdata

head(counts(d))
samples(d)
head(names(d))
length(d)
d[1:20, ]
d[1:20, 1:3]

### Filtering
# Check what is the minimal number of replicates per condition 
table(samples(d)$group)
d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3, 
  min_samps_feature_prop = 0)
plotData(d)

### Calculate dispersion
d <- dmDispersion(d, verbose = 0, BPPARAM = BiocParallel::SerialParam())
plotDispersion(d)

head(mean_expression(d))
common_dispersion(d)
head(genewise_dispersion(d))

### Fit full model proportions
d <- dmFit(d, verbose = 0, BPPARAM = BiocParallel::SerialParam())

head(proportions(d))
head(statistics(d))

### Fit null model proportions and test for DS
d <- dmTest(d, BPPARAM = BiocParallel::SerialParam())
plotTest(d)

head(proportions(d))
head(statistics(d))
head(results(d))

### Plot feature proportions for top DS gene
res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

gene_id <- res$gene_id[1]

plotFit(d, gene_id = gene_id)
plotFit(d, gene_id = gene_id, plot_type = "lineplot")
plotFit(d, gene_id = gene_id, plot_type = "ribbonplot")












