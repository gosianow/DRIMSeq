context("DM adjustement")


gamma0 = 10
y = matrix(c(30, 75, 35, 70), nrow = 2)
pi = c(0.3, 0.7)

test_that("dm_adjustmentOneGeneOneGroup returns the right values", {
  
  expect_equal(round(dm_adjustmentOneGeneOneGroup(y, gamma0, pi), 4), 2.6562)
  
})




