context("DM likelihood and score")

gamma0 = 10
y = matrix(c(35, 70, 100, 100), nrow = 2)
pi = 0.3


test_that("dm_lik and dm_likG are equal", {
  expect_equal(dm_lik(pi, gamma0, y), dm_likG(pi, gamma0, y))
})


test_that("dm_score and dm_scoreG are equal", {
  expect_equal(dm_score(pi, gamma0, y), dm_scoreG(pi, gamma0, y))
})


y = matrix(c(0, 0, 35, 70), nrow = 2)

test_that("dm_lik and dm_likG are diff. for matrix with 0 column", {
  expect_true(dm_lik(pi, gamma0, y) < dm_likG(pi, gamma0, y))
})

