test_that("show", {
  data(ChainSC)
  data(ChainSCReg)
  expect_output(show(ChainSC), "An object of class BASiCS_Chain")
  expect_output(show(ChainSCReg), "An object of class BASiCS_Chain")
  expect_output(show(Summary(ChainSC)), "An object of class BASiCS_Summary")
  d <- BASiCS_TestDE(ChainSCReg, ChainSCReg)
  expect_output(show(d), "An object of class BASiCS_ResultsDE, containing")
  v <- BASiCS_DetectHVG(ChainSCReg, PercentileThreshold = 0.9)
  expect_output(show(v), "An object of class BASiCS_ResultVG.")
})

test_that("Summary na.rm", {
  data(ChainSC)
  ChainSC@parameters$nu[1, 1] <- NA
  ChainSC@parameters$nu[1, 4] <- NA
  s1 <- Summary(ChainSC, na.rm = TRUE)
  s2 <- Summary(ChainSC, na.rm = FALSE)
  expect_false(any(is.na(s1@parameters$nu)))
  expect_true(all(is.na(s2@parameters$nu[1, ])))
  expect_true(all(is.na(s2@parameters$nu[4, ])))
})

test_that("subset", {
  data(ChainSC)
  data(ChainSCReg)
  sc <- subset(
    ChainSC,
    Genes = rownames(ChainSC)[1:3],
    Cells = colnames(ChainSC)[1:2],
    Iterations = 1:10
  )
  expect_equal(dim(sc@parameters$mu), c(10, 3))
  expect_equal(dim(sc@parameters$nu), c(10, 2))
  sc <- subset(
    ChainSCReg,
    Genes = rownames(ChainSC)[1:3],
    Cells = colnames(ChainSC)[1:2],
    Iterations = 1:10
  )
  expect_equal(dim(sc@parameters$epsilon), c(10, 3))
  expect_equal(dim(sc@parameters$nu), c(10, 2))
})

test_that("subset respects order", {
  data(ChainSC)
  data(ChainSCReg)
  sc <- subset(
    ChainSC,
    Genes = rownames(ChainSC)[3:1]
  )
  expect_equal(rownames(sc), rownames(ChainSC)[3:1])
  sc <- subset(
    ChainSC,
    Genes = 3:1
  )
  expect_equal(rownames(sc), rownames(ChainSC)[3:1])
  ind <- rep(c(TRUE, FALSE), length.out = nrow(ChainSC))
  sc <- subset(
    ChainSC,
    Genes = ind
  )
  expect_equal(rownames(sc), rownames(ChainSC)[ind])
})



test_that("subset with reffreq", {
  Data <- makeExampleBASiCS_Data(WithBatch = TRUE)
  Chain <- run_MCMC(
    Data, N = 50, Thin = 2, Burn = 10,   Regression = TRUE,
    PrintProgress = FALSE, WithSpikes = FALSE
  )
  expect_error(subset(Chain, Genes = c("Gene1", "Gene2")), NA)
})

test_that("dimnames", {
  data(ChainSC)
  expect_equal(colnames(ChainSC), colnames(ChainSC@parameters$nu))
  expect_equal(rownames(ChainSC), colnames(ChainSC@parameters$mu))
  expect_equal(dimnames(ChainSC), list(rownames(ChainSC), colnames(ChainSC)))
})

test_that("subsets for BASiCS_ResultsDE", {
  data(ChainSC)
  data(ChainSCReg)
  d <- BASiCS_TestDE(ChainSCReg, ChainSCReg)
  g <- d[[1]]@Table$GeneName[1:10]
  ds <- d[g, ]
  expect_equal(ds[[1]]@Table, d[[1]]@Table[1:10, ])
})
