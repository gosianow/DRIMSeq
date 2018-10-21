context("total switch")

test_that("total switch retains precision estimates", {

  cts <- matrix(rnbinom(100, mu=100, size=10), ncol=10,
                dimnames=list(1:10,paste0("s",1:10)))
  cts[1,1:5] <- 0
  cts[2,6:10] <- 0

  counts <- data.frame(gene_id=factor(rep(1:5, each=2)),
                       feature_id=factor(1:10), cts)
  samples <- data.frame(sample_id=paste0("s",1:10), condition=rep(c("A","B"),each=5))
  d <- dmDSdata(counts, samples)
  design <- model.matrix(~condition, samples(d))

  d <- dmPrecision(d, design=design, prec_subset=1)
  expect_true(is.na(genewise_precision(d)$genewise_precision[1]))

  d2 <- dmPrecision(d, design=design, prec_subset=1, add_uniform=TRUE)
  expect_true(!is.na(genewise_precision(d2)$genewise_precision)[1])
  
  d <- dmFit(d, design=design)
  d <- dmTest(d, coef="conditionB")
  res <- results(d)
  expect_true(is.na(res$pvalue[1]))

  d2 <- dmFit(d2, design=design, add_uniform=TRUE)
  d2 <- dmTest(d2, coef="conditionB")
  res <- results(d2)
  expect_true(!is.na(res$pvalue[1]))

})
