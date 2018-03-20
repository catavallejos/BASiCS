context("Differential test\n")

test_that("Differential test is correct", 
{
  data(ChainSC)
  data(ChainRNA)
            
  Test <- BASiCS_TestDE(Chain1 = ChainSC, Chain2 = ChainRNA,
                        GroupLabel1 = 'SC', GroupLabel2 = 'P&S',
                        EpsilonM = log2(1.5), EpsilonD = log2(1.5), 
                        OffSet = TRUE, Plot = FALSE, PlotOffset = FALSE)
            
  # Classification frequency
            
  FreqMean0 <- c(369,  22,   9)
  FreqMean <- as.vector(table(Test$TableMean$ResultDiffMean))
  expect_that(all.equal(FreqMean, FreqMean0), is_true())
            
  FreqDisp0 <- c(296,   2,  71)
  FreqDisp <- as.vector(table(Test$TableDisp$ResultDiffDisp))
  expect_that(all.equal(FreqDisp, FreqDisp0), is_true())
            
  # Posterior probabilities
            
  ProbMean0 <- c(0.00, 0.01, 0.87, 0.48, 0.28, 
                 0.00, 0.07, 0.28, 0.09, 0.01)
  ProbMean <- round(Test$TableMean$ProbDiffMean[1:10], 2)
  expect_that(all.equal(ProbMean, ProbMean0), is_true())
            
  ProbDisp0 <- c(0.83, 0.69, 0.93, 0.75, 0.71, 
                 0.59, 0.56, 0.88, 0.60, 0.64)
  ProbDisp <- round(tail(Test$TableDisp$ProbDiffDisp, 10), 2)
  expect_that(all.equal(ProbDisp, ProbDisp0), is_true())
            
  # Log2 fold changes
            
  Lfc2Mean0 <- c(-0.22, -0.08, -0.84, -0.57, -0.38,  
                  0.01, -0.02,  0.48, -0.32,  0.13)
  Lfc2Mean <- round(Test$TableMean$MeanLog2FC[1:10], 2)
  expect_that(all.equal(Lfc2Mean, Lfc2Mean0), is_true())
            
  Lfc2Disp0 <- c( 1.18, -0.82,  2.17,  0.96, -0.98,  
                  0.72,  0.59,  1.74,  0.52,  0.69)
  Lfc2Disp <- round(tail(Test$TableDisp$DispLog2FC, 10),2)
  expect_that(all.equal(Lfc2Disp, Lfc2Disp0), is_true())
            
})




