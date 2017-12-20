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
            
  FreqMean0 <- c(367,  17,  16)
  FreqMean <- as.vector(table(Test$TableMean$ResultDiffMean))
  expect_that(all.equal(FreqMean, FreqMean0), is_true())
            
  FreqDisp0 <- c(293,  74)
  FreqDisp <- as.vector(table(Test$TableDisp$ResultDiffDisp))
  expect_that(all.equal(FreqDisp, FreqDisp0), is_true())
            
  # Posterior probabilities
            
  ProbMean0 <- c(1.00, 1.00, 1.00, 1.00, 1.00, 
                 1.00, 1.00, 1.00, 1.00 ,0.97)
  ProbMean <- round(Test$TableMean$ProbDiffMean[1:10], 2)
  expect_that(all.equal(ProbMean, ProbMean0), is_true())
            
  ProbDisp0 <- c(0.41, 0.40, 0.40, 0.39, 0.39, 
                 0.37, 0.36, 0.36, 0.32, 0.31)
  ProbDisp <- round(tail(Test$TableDisp$ProbDiffDisp, 10), 2)
  expect_that(all.equal(ProbDisp, ProbDisp0), is_true())
            
  # Log2 fold changes
            
  Lfc2Mean0 <- c(-1.53,  1.71, -4.07, -4.27,  1.95, 
                 -1.00,  1.29,  0.90,  1.65,  1.43)
  Lfc2Mean <- round(Test$TableMean$MeanLog2FC[1:10], 2)
  expect_that(all.equal(Lfc2Mean, Lfc2Mean0), is_true())
            
  Lfc2Disp0 <- c( 0.23,  0.12, -0.16, -0.16,  0.34, 
                 -0.30, -0.06,  0.12,  0.24,  0.23)
  Lfc2Disp <- round(tail(Test$TableDisp$DispLog2FC, 10),2)
  expect_that(all.equal(Lfc2Disp, Lfc2Disp0), is_true())
            
})




