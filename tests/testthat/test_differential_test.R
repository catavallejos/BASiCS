context("Differential test\n")

test_that("Differential test is correct", {
  data(ChainSC)
  data(ChainRNA)
            
  Test <- BASiCS_TestDE(Chain1 = ChainSC, Chain2 = ChainRNA,
                        GroupLabel1 = 'SC', GroupLabel2 = 'P&S',
                        EpsilonM = log2(1.5), EpsilonD = log2(1.5), 
                        OffSet = TRUE, Plot = FALSE)
            
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
