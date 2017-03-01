context("DM likelihood and score")

prec = 10
y = matrix(c(35, 70, 100, 100), nrow = 2)
prop = 0.3


test_that("dm_lik and dm_likG are equal", {
  expect_equal(dm_lik(prop, prec, y), dm_likG(prop, prec, y))
})


test_that("dm_score and dm_scoreG are equal", {
  expect_equal(dm_score(prop, prec, y), dm_scoreG(prop, prec, y))
})


y = matrix(c(0, 0, 35, 70), nrow = 2)

test_that("dm_lik and dm_likG are diff. for matrix with 0 column", {
  expect_true(dm_lik(prop, prec, y) < dm_likG(prop, prec, y))
})

