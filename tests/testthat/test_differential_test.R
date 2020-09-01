context("Differential test")

data(ChainSC)
data(ChainRNA)

test_that("Differential test is correct", {

  Test <- BASiCS_TestDE(
    Chain1 = ChainSC,
    Chain2 = ChainRNA,
    GroupLabel1 = 'SC',
    GroupLabel2 = 'P&S',
    EpsilonM = log2(1.5),
    EpsilonD = log2(1.5),
    MinESS = NA,
    OffSet = TRUE,
    Plot = FALSE
  )
            
  # Classification frequency
            
  FreqMean0 <- c(335,  8,   7)
  FreqMean <- as.vector(table(Test@Results$Mean@Table$ResultDiffMean))
  expect_equal(FreqMean, FreqMean0)
            
  FreqDisp0 <- c(15, 295,   1,  39)
  FreqDisp <- as.vector(table(Test@Results$Disp@Table$ResultDiffDisp))
  expect_equal(FreqDisp, FreqDisp0)
            
  # Posterior probabilities  
  ProbMean0 <- c(0.00, 0.01, 0.77, 0.43, 0.19)
  ProbMean <- round(Test@Results$Mean@Table$ProbDiffMean[1:5], 2)
  expect_equal(ProbMean, ProbMean0)
            
  ProbDisp0 <- c(0.69, 0.60, 0.60, 1.00, 0.40)
  ProbDisp <- round(tail(Test@Results$Disp@Table$ProbDiffDisp, 5), 2)
  expect_equal(ProbDisp, ProbDisp0)
            
  # Log2 fold changes      
  Lfc2Mean0 <- c(-0.17, -0.02, -0.78, -0.51, -0.32)
  Lfc2Mean <- round(Test@Results$Mean@Table$MeanLog2FC[1:5], 2)
  expect_equal(Lfc2Mean, Lfc2Mean0)
            
  Lfc2Disp0 <- c(0.62,  0.67,  0.54, 4.39, -0.12)
  Lfc2Disp <- round(tail(Test@Results$Disp@Table$DispLog2FC, 5),2)
  expect_equal(Lfc2Disp, Lfc2Disp0)
})


test_that("Differential test fails when different number of samples", {
  Chain2 <- ChainSC
  Chain2@parameters <- lapply(Chain2@parameters, function(x) x[1:10, ])
  expect_error(
    BASiCS_TestDE(ChainSC, Chain2),
    "Chains must have an equal number of samples to run BASiCS_TestDE."
  )
})

test_that("GeneSelect functions reasonably well", {

  # Test usage of GenesSelect
  data(ChainSC)
  data(ChainRNA)

  nGenes <- ncol(ChainSC@parameters$mu)
  set.seed(10)
  GenesSelect <- rbinom(n = nGenes, size = 1, prob = 0.9)
  GenesSelect <- as.logical(GenesSelect)
  Test <- BASiCS_TestDE(
    Chain1 = ChainSC,
    Chain2 = ChainRNA,
    GroupLabel1 = 'SC',
    GroupLabel2 = 'P&S',
    GenesSelect = GenesSelect,
    EpsilonM = log2(1.5),
    EpsilonD = log2(1.5),
    MinESS = NA,
    OffSet = TRUE,
    Plot = FALSE,
    PlotOffset = FALSE
  )
            
  # Classification frequency
  FreqMean0 <- c(26, 311,   6,   7)
  FreqMean <- as.vector(table(Test@Results$Mean@Table$ResultDiffMean))
  expect_equal(FreqMean, FreqMean0)
            
  FreqDisp0 <- c(26,  13, 279,   1,  31)
  FreqDisp <- as.vector(table(Test@Results$Disp@Table$ResultDiffDisp))
  expect_equal(FreqDisp, FreqDisp0)
            
  # Posterior probabilities
  ProbMean0 <- c(0.00, 0.01, 0.77, 0.43, 0.19)
  ProbMean <- round(Test@Results$Mean@Table$ProbDiffMean[1:5], 2)
  expect_equal(ProbMean, ProbMean0)
            
  ProbDisp0 <- c(0.69, 0.60, 0.60, 1.00, 0.40)
  ProbDisp <- round(tail(Test@Results$Disp@Table$ProbDiffDisp, 5), 2)
  expect_equal(ProbDisp, ProbDisp0)
            
  # Log2 fold changes
  Lfc2Mean0 <- c(-0.17, -0.02, -0.78, -0.51, -0.32)
  Lfc2Mean <- round(Test@Results$Mean@Table$MeanLog2FC[1:5], 2)
  expect_equal(Lfc2Mean, Lfc2Mean0)
            
  Lfc2Disp0 <- c(0.62,  0.67,  0.54, 4.39, -0.12)
  Lfc2Disp <- round(tail(Test@Results$Disp@Table$DispLog2FC, 5),2)
  expect_equal(Lfc2Disp, Lfc2Disp0)
})


test_that("CheckESS works", {
  data(ChainSC)
  data(ChainRNA)

  MinESS <- 100
  Test <- BASiCS_TestDE(
    Chain1 = ChainSC,
    Chain2 = ChainRNA,
    EpsilonM = log2(1.5),
    EpsilonD = log2(1.5),
    OffSet = TRUE,
    Plot = FALSE,
    MinESS = MinESS,
    PlotOffset = FALSE
  )
  
  me1 <- coda::effectiveSize(ChainSC@parameters$mu)
  me2 <- coda::effectiveSize(ChainRNA@parameters$mu)
  # Classification frequency
  expect_equivalent(
    Test@Results$Mean@Table$ResultDiffMean == "ExcludedLowESS",
    !(me1 > MinESS & me2 > MinESS)
  )

  md1 <- coda::effectiveSize(ChainSC@parameters$delta)
  md2 <- coda::effectiveSize(ChainRNA@parameters$delta)
  # Classification frequency
  expect_equivalent(
    Test@Results$Disp@Table$ResultDiffDisp == "ExcludedLowESS",
    !(md1 > MinESS & md2 > MinESS)
  )
})



test_that("EpsilonM = 0 case (reg)", {
  data(ChainSC)
  data(ChainRNA)

  Test <- BASiCS_TestDE(
    Chain1 = ChainSC,
    Chain2 = ChainRNA,
    GroupLabel1 = "SC",
    GroupLabel2 = "P&S",
    EpsilonM = 0,
    EpsilonD = 0,
    OffSet = TRUE,
    Plot = FALSE,
    MinESS = NA,
    PlotOffset = FALSE
  )
  expect_equal(
    as.numeric(table(Test@Results$Mean@Table$ResultDiffMean)),
    c(180, 98, 72)
  )
  expect_equal(
    as.numeric(table(Test@Results$Disp@Table$ResultDiffDisp)),
    c(170, 147, 33)
  )
})


test_that("OffSet=FALSE works", {
  data(ChainSC)
  data(ChainRNA)
  expect_error(
    BASiCS_TestDE(
      Chain1 = ChainSC,
      Chain2 = ChainRNA,
      OffSet = FALSE
    ),
    NA
  )
})

test_that("Utility methods work", {
  data(ChainSC)
  data(ChainRNA)

  Test <- BASiCS_TestDE(
    Chain1 = ChainSC,
    Chain2 = ChainRNA,
    GroupLabel1 = "SC",
    GroupLabel2 = "P&S",
    EpsilonM = 0,
    EpsilonD = 0,
    OffSet = TRUE,
    Plot = FALSE,
    MinESS = NA,
    PlotOffset = FALSE
  )
  expect_is(Test, "BASiCS_ResultsDE")
  expect_is(Test[["Mean"]], "BASiCS_ResultDE")
  f <- Test[["Mean"]][1:10, ]
  expect_warning(df <- format(f, Filter = FALSE))
  expect_equal(
    Test@Results$Mean@Table[1:10, "GeneName"],
    as.character(df[["GeneName"]])
  )
  d <- rowData(Test)
  ind <- seq_len(nrow(d))
  d$ind <- ind
  rowData(Test) <- d
  expect_equal(rowData(Test)$ind, ind)
  df <- as.data.frame(Test, Parameter = "Mean")
  expect_is(df, "data.frame")
  p <- 0.9
  df <- as.data.frame(Test, Parameter = "Mean", ProbThreshold = p)
  expect_true(all(df$ProbDiffMean > p))
})
