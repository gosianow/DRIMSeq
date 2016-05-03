context("DM fit proportions")

################################################################
# dm_fitOneGeneOneGroup
################################################################

gamma0 = 10
y = matrix(c(30, 70, 60, 140), nrow = 2)


test_that("dm_fitOneGeneOneGroup returns list(pi = pi, stats = c(lik = lik, df = df))", {
  
  dm_fitOneGeneOneGroup_out <- dm_fitOneGeneOneGroup(y, gamma0, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE)
  
  expect_is(dm_fitOneGeneOneGroup_out, "list")
  expect_true(all(names(dm_fitOneGeneOneGroup_out) == c("pi", "stats")))
  expect_true(all(names(dm_fitOneGeneOneGroup_out$stats) == c("lik", "df")))
  
})


gamma0 = 10
y = matrix(c(30, 70, 0, 0), nrow = 2)


test_that("dm_fitOneGeneOneGroup returns rowSums(y)/sum(y) for y with one non-zero column", {
  
  expect_equal(dm_fitOneGeneOneGroup(y, gamma0, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE)$pi, c(0.3, 0.7))
  
})



gamma0 = 10
y = matrix(c(30, 0, 30, 0), nrow = 2)


test_that("dm_fitOneGeneOneGroup returns NAs for y with one non-zero row", {
  
  expect_equal(dm_fitOneGeneOneGroup(y, gamma0, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE)$pi, c(NA, NA))
  
})



################################################################
# dm_fitOneGeneManyGroups
################################################################



gamma0 = 10
y = matrix(c(30, 70, 60, 140, 50, 50, 100, 100), nrow = 2)
ngroups = 2
lgroups = c("C1", "C2")
igroups = list(C1 = 1:2, C2 = 3:4)



test_that("dm_fitOneGeneOneGroup returns list(pi = pi, stats = stats)", {
  
  dm_fitOneGeneManyGroups_out <- dm_fitOneGeneManyGroups(y, ngroups, lgroups, igroups, gamma0, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE)
  
  expect_is(dm_fitOneGeneManyGroups_out, "list")
  expect_true(all(names(dm_fitOneGeneManyGroups_out) == c("pi", "stats")))
  expect_true(all(names(dm_fitOneGeneManyGroups_out$stats) == c("C1", "C2")))
  expect_true(all(colnames(dm_fitOneGeneManyGroups_out$pi) == c("C1", "C2")))
  
  
})




gamma0 = 10
y = matrix(c(30, 70, 60, 140, 0, 50, 0, 100), nrow = 2)
ngroups = 2
lgroups = c("C1", "C2")
igroups = list(C1 = 1:2, C2 = 3:4)



test_that("dm_fitOneGeneOneGroup returns NAs when any of the groups has NAs", {
  
  dm_fitOneGeneManyGroups_out <- dm_fitOneGeneManyGroups(y, ngroups, lgroups, igroups, gamma0, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE)
  
  expect_true(all(is.na(dm_fitOneGeneManyGroups_out$pi)))
  expect_true(all(is.na(dm_fitOneGeneManyGroups_out$stats)))
  
})


