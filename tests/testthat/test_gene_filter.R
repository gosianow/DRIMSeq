context("gene filter")

test_that("gene filtering is working", {

  cts <- matrix(rnbinom(120, mu=100, size=10), ncol=10,
                dimnames=list(1:12,paste0("s",1:10)))
  cts[1,] <- 0
  cts[2:3,] <- 10
  cts[1:3,1] <- c(10,0,0)
  counts <- data.frame(gene_id=factor(rep(1:4, each=3)),
                       feature_id=factor(1:12), cts)
  samples <- data.frame(sample_id=paste0("s",1:10), condition=rep(c("A","B"),each=5))
  d <- dmDSdata(counts, samples)
  d <- dmFilter(d, min_samps_gene_expr = 10, min_gene_expr = 10,
                min_samps_feature_expr = 5, min_feature_expr = 0,
                min_samps_feature_prop = 5, min_feature_prop = 0)
  
  d2 <- dmFilter(d, min_samps_gene_expr = 10, min_gene_expr = 10,
                 min_samps_feature_expr = 5, min_feature_expr = 10,
                 min_samps_feature_prop = 5, min_feature_prop = 0)
  expect_true(length(d2) == 4)
  sum(counts(d2[1,])$s1) # has 0 counts for this gene

  # now use 'run_gene_twice'
  d3 <- dmFilter(d, min_samps_gene_expr = 10, min_gene_expr = 10,
                 min_samps_feature_expr = 5, min_feature_expr = 10,
                 min_samps_feature_prop = 5, min_feature_prop = 0,
                 run_gene_twice=TRUE)
  expect_true(length(d3) == 3)

})
