context("MatrixList")

x1 <- matrix(1, 10, 6)
colnames(x1) <- paste0("C", 1:6)
rownames(x1) <- paste0("R1", 1:10)
x2 <- matrix(2, 5, 6)
colnames(x2) <- paste0("C", 1:6)
rownames(x2) <- paste0("R2", 1:5)

x <- MatrixList(x1 = x1, x2 = x2)


test_that("methods return correct arributes of MatrixList", {
  expect_equal(names(x), c("x1", "x2"))
  expect_equal(rownames(x), c(paste0("R1", 1:10), paste0("R2", 1:5)))
  expect_equal(colnames(x), paste0("C", 1:6))
  expect_equal(length(x), 2)
  expect_equal(elementLengths(x), c(x1 = 10, x2 = 5))
  expect_equal(dim(x), c(15, 6))
  expect_equal(nrow(x), 15)
  expect_equal(ncol(x), 6)
})


test_that("subsetting of MatrixList is correct", {
  expect_equal(x[["x2"]], x2)
  expect_equal(x$x2, x2)
  expect_equal(x[1, 1], MatrixList(x1 = x1[, 1, drop = FALSE]))
})



