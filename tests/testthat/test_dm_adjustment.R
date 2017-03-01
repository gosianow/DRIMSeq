context("DM adjustement")


prec = 10
y = matrix(c(30, 75, 35, 70), nrow = 2)
prop = c(0.3, 0.7)

test_that("dm_CRadjustmentOneGroup returns the right values", {
  
  expect_equal(round(dm_CRadjustmentOneGroup(y, prec, prop), 4), 2.6562)
  
})




