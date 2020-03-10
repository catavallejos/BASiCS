context("Differential test (regression case)")

test_that("Differential test is correct (regression case)", {
  data(ChainSCReg)
  data(ChainRNAReg)

  Test <- BASiCS_TestDE(
    Chain1 = ChainSCReg,
    Chain2 = ChainRNAReg,
    GroupLabel1 = "SC",
    GroupLabel2 = "P&S",
    EpsilonM = log2(1.5),
    EpsilonD = log2(1.5),
    EpsilonR = log2(1.5) / log2(exp(1)),
    MinESS = NA,
    OffSet = TRUE, 
    Plot = FALSE,
    PlotOffset = FALSE
  )

  # Classification frequency

  FreqMean0 <- c(333, 8, 9)
  FreqMean <- as.vector(table(Test$TableMean$ResultDiffMean))
  expect_equal(FreqMean, FreqMean0)

  FreqDisp0 <- c(17, 233, 100)
  FreqDisp <- as.vector(table(Test$TableDisp$ResultDiffDisp))
  expect_equal(FreqDisp, FreqDisp0)

  FreqRes0 <- c(338, 4, 8)
  FreqRes <- as.vector(table(Test$TableResDisp$ResultDiffResDisp))
  expect_equal(FreqRes, FreqRes0)

  # Posterior probabilities

  ProbMean0 <- c(0.01, 0.01, 0.75, 0.32, 0.28)
  ProbMean <- round(Test$TableMean$ProbDiffMean[1:5], 2)
  expect_equal(ProbMean, ProbMean0)

  ProbDisp0 <- c(0.72, 0.72, 0.72, 1.00, 0.59)
  ProbDisp <- round(tail(Test$TableDisp$ProbDiffDisp, 5), 2)
  expect_equal(ProbDisp, ProbDisp0)

  ProbRes0 <- c(0.60, 0.56, 0.63, 1.00, 0.52)
  ProbRes <- round(tail(Test$TableRes$ProbDiffRes, 5), 2)
  expect_equal(ProbRes, ProbRes0)

  # Log2 fold changes
  Lfc2Mean0 <- c(-0.17, -0.02, -0.86, -0.48, -0.40)
  Lfc2Mean <- round(Test$TableMean$MeanLog2FC[1:5], 2)
  expect_equal(Lfc2Mean, Lfc2Mean0)

  Lfc2Disp0 <- c(0.92, 1.02, 0.86, 6.23, 0.26)
  Lfc2Disp <- round(tail(Test$TableDisp$DispLog2FC, 5), 2)
  expect_equal(Lfc2Disp, Lfc2Disp0)
})

test_that("Differential test requires same regression setting", {
  data(ChainSCReg)
  data(ChainRNA)
  expect_error(
    BASiCS_TestDE(
      Chain1 = ChainSCReg,
      Chain2 = ChainRNA,
      Plot = FALSE,
      PlotOffset = FALSE
    ),
    "Both chains should be run with the same setting for Regression."
  )
})

test_that("Differential test when testing ESS", {
  data(ChainSCReg)
  data(ChainRNA)
  expect_error(
    BASiCS_TestDE(
      Chain1 = ChainSCReg,
      Chain2 = ChainRNAReg,
      Plot = FALSE,
      PlotOffset = FALSE
    ),
    NA
  )
  })

test_that("CheckESS works", {
  data(ChainSCReg)
  data(ChainRNAReg)

  MinESS <- 100
  Test <- BASiCS_TestDE(
    Chain1 = ChainSCReg,
    Chain2 = ChainRNAReg,
    EpsilonM = log2(1.5),
    EpsilonD = log2(1.5),
    OffSet = TRUE,
    Plot = FALSE,
    MinESS = MinESS,
    PlotOffset = FALSE
  )
  ee1 <- coda::effectiveSize(ChainSCReg@parameters$epsilon)
  ee2 <- coda::effectiveSize(ChainRNAReg@parameters$epsilon)

  # Classification frequency
  expect_equivalent(
    Test$TableResDisp$ResultDiffResDisp == "ExcludedLowESS",
    !(ee1 > MinESS & ee2 > MinESS)
  )
})



test_that("EpsilonM = 0 case (reg)", {
  data(ChainSCReg)
  data(ChainRNAReg)

  Test <- BASiCS_TestDE(
    Chain1 = ChainSCReg,
    Chain2 = ChainRNAReg,
    GroupLabel1 = "SC",
    GroupLabel2 = "P&S",
    EpsilonM = 0,
    EpsilonD = 0,
    EpsilonR = 0,
    OffSet = TRUE,
    Plot = FALSE,
    MinESS = NA,
    PlotOffset = FALSE
  )
  expect_equal(
    as.numeric(table(Test$TableMean$ResultDiffMean)),
    c(150, 119, 81)
  )
  expect_equal(
    as.numeric(table(Test$TableDisp$ResultDiffDisp)),
    c(200, 98, 52)
  )
  expect_equal(
    as.numeric(table(Test$TableResDisp$ResultDiffResDisp)),
    c(333, 6, 11)
  )
})
