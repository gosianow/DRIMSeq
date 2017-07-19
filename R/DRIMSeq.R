#' @import BiocGenerics
#' @import methods
#' @import BiocParallel
#' @import edgeR
#' @import GenomicRanges
NULL

# ### EXAMPLES ###
# # --------------------------------------------------------------------------
# # Create dmDSdata object 
# # --------------------------------------------------------------------------
# ## Get kallisto transcript counts from the 'PasillaTranscriptExpr' package
# 
# library(PasillaTranscriptExpr)
# 
# data_dir  <- system.file("extdata", package = "PasillaTranscriptExpr")
# 
# ## Load metadata
# pasilla_metadata <- read.table(file.path(data_dir, "metadata.txt"), 
#   header = TRUE, as.is = TRUE)
# 
# ## Load counts
# pasilla_counts <- read.table(file.path(data_dir, "counts.txt"), 
#   header = TRUE, as.is = TRUE)
# 
# ## Create a pasilla_samples data frame
# pasilla_samples <- data.frame(sample_id = pasilla_metadata$SampleName, 
#   group = pasilla_metadata$condition)
# levels(pasilla_samples$group)
# 
# ## Create a dmDSdata object
# d <- dmDSdata(counts = pasilla_counts, samples = pasilla_samples)
# 
# ## Use a subset of genes, which is defined in the following file
# gene_id_subset <- readLines(file.path(data_dir, "gene_id_subset.txt"))
# 
# d <- d[names(d) %in% gene_id_subset, ]
# 
# # --------------------------------------------------------------------------
# # Differential transcript usage analysis - simple two group comparison 
# # --------------------------------------------------------------------------
# 
# ## Filtering
# ## Check what is the minimal number of replicates per condition 
# table(samples(d)$group)
# 
# d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3,
#   min_gene_expr = 10, min_feature_expr = 10)
# 
# plotData(d)
# 
# ## Create the design matrix
# design_full <- model.matrix(~ group, data = samples(d))
# 
# ## To make the analysis reproducible
# set.seed(123)
# ## Calculate precision
# d <- dmPrecision(d, design = design_full)
# 
# plotPrecision(d)
# 
# head(mean_expression(d))
# common_precision(d)
# head(genewise_precision(d))
# 
# ## Fit full model proportions
# d <- dmFit(d, design = design_full)
# 
# ## Get fitted proportions
# head(proportions(d))
# ## Get the DM regression coefficients (gene-level) 
# head(coefficients(d))
# ## Get the BB regression coefficients (feature-level) 
# head(coefficients(d), level = "feature")
# 
# ## Fit null model proportions and perform the LR test to detect DTU
# d <- dmTest(d, coef = "groupKD")
# 
# ## Plot the gene-level p-values
# plotPValues(d)
# 
# ## Get the gene-level results
# head(results(d))
# 
# ## Plot feature proportions for a top DTU gene
# res <- results(d)
# res <- res[order(res$pvalue, decreasing = FALSE), ]
# 
# top_gene_id <- res$gene_id[1]
# 
# plotProportions(d, gene_id = top_gene_id, group_variable = "group")
# 
# plotProportions(d, gene_id = top_gene_id, group_variable = "group", 
#   plot_type = "lineplot")
# 
# plotProportions(d, gene_id = top_gene_id, group_variable = "group", 
#   plot_type = "ribbonplot")
# 
# 
# 
# 
# 
# 
# # --------------------------------------------------------------------------
# # Create dmSQTLdata object
# # --------------------------------------------------------------------------
# # Use subsets of data defined in the GeuvadisTranscriptExpr package
# 
# library(GeuvadisTranscriptExpr)
# 
# geuv_counts <- GeuvadisTranscriptExpr::counts
# geuv_genotypes <- GeuvadisTranscriptExpr::genotypes
# geuv_gene_ranges <- GeuvadisTranscriptExpr::gene_ranges
# geuv_snp_ranges <- GeuvadisTranscriptExpr::snp_ranges
# 
# colnames(geuv_counts)[c(1,2)] <- c("feature_id", "gene_id")
# colnames(geuv_genotypes)[4] <- "snp_id"
# geuv_samples <- data.frame(sample_id = colnames(geuv_counts)[-c(1,2)])
# 
# d <- dmSQTLdata(counts = geuv_counts, gene_ranges = geuv_gene_ranges,
#   genotypes = geuv_genotypes, snp_ranges = geuv_snp_ranges, 
#   samples = geuv_samples, window = 5e3)
# 
# # --------------------------------------------------------------------------
# # sQTL analysis - simple group comparison
# # --------------------------------------------------------------------------
# 
# ## Filtering
# d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_expr = 5,
#   minor_allele_freq = 5, min_gene_expr = 10, min_feature_expr = 10)
# 
# plotData(d)
# 
# ## To make the analysis reproducible
# set.seed(123)
# ## Calculate precision
# d <- dmPrecision(d)
# 
# plotPrecision(d)
# 
# ## Fit full model proportions
# d <- dmFit(d)
# 
# ## Fit null model proportions, perform the LR test to detect tuQTLs 
# ## and use the permutation approach to adjust the p-values
# d <- dmTest(d)
# 
# ## Plot the gene-level p-values
# plotPValues(d)
# 
# ## Get the gene-level results
# head(results(d))



