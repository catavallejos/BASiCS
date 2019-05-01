context("Differential test\n")

data(ChainSC)
data(ChainRNA)

test_that("Differential test is correct", {
            
  Test <- BASiCS_TestDE(Chain1 = ChainSC, Chain2 = ChainRNA,
                        GroupLabel1 = 'SC', GroupLabel2 = 'P&S',
                        EpsilonM = log2(1.5), EpsilonD = log2(1.5), 
                        OffSet = TRUE, Plot = FALSE, PlotOffset = FALSE)
            
  # Classification frequency
            
  FreqMean0 <- c(335,  10,   5)
  FreqMean <- as.vector(table(Test$TableMean$ResultDiffMean))
  expect_that(all.equal(FreqMean, FreqMean0), is_true())
            
  FreqDisp0 <- c(15, 295,   1,  39)
  FreqDisp <- as.vector(table(Test$TableDisp$ResultDiffDisp))
  expect_that(all.equal(FreqDisp, FreqDisp0), is_true())
            
  # Posterior probabilities
            
  ProbMean0 <- c(0.00, 0.01, 0.87, 0.49, 0.28)
  ProbMean <- round(Test$TableMean$ProbDiffMean[1:5], 2)
  expect_that(all.equal(ProbMean, ProbMean0), is_true())
            
  ProbDisp0 <- c(0.69, 0.60, 0.60, 1.00, 0.40)
  ProbDisp <- round(tail(Test$TableDisp$ProbDiffDisp, 5), 2)
  expect_that(all.equal(ProbDisp, ProbDisp0), is_true())
            
  # Log2 fold changes
            
  Lfc2Mean0 <- c(-0.23, -0.08, -0.84, -0.57, -0.38)
  Lfc2Mean <- round(Test$TableMean$MeanLog2FC[1:5], 2)
  expect_that(all.equal(Lfc2Mean, Lfc2Mean0), is_true())
            
  Lfc2Disp0 <- c(0.62,  0.67,  0.54, 4.39, -0.12)
  Lfc2Disp <- round(tail(Test$TableDisp$DispLog2FC, 5),2)
  expect_that(all.equal(Lfc2Disp, Lfc2Disp0), is_true())
            
})

test_that("Differential test fails when different number of samples", {
  Chain2 <- ChainSC
  Chain2@parameters <- lapply(Chain2@parameters, function(x) x[1:10, ])
  expect_error(
    BASiCS_TestDE(ChainSC, Chain2),
    "Chains must have an equal number of samples to run BASiCS_TestDE."
  )
})


