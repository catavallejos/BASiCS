context("Differential test (regression case)\n")

test_that("Differential test is correct (regression case)", 
{
  data(ChainSCReg)
  data(ChainRNAReg)
            
  Test <- BASiCS_TestDE(Chain1 = ChainSCReg, Chain2 = ChainRNAReg,
                        GroupLabel1 = 'SC', GroupLabel2 = 'P&S',
                        EpsilonM = log2(1.5), EpsilonD = log2(1.5), 
                        EpsilonR = log2(1.5)/log2(exp(1)),
                        OffSet = TRUE, Plot = FALSE, PlotOffset = FALSE)
            
  # Classification frequency
            
  FreqMean0 <- c(324,  18,   8)
  FreqMean <- as.vector(table(Test$TableMean$ResultDiffMean))
  expect_that(all.equal(FreqMean, FreqMean0), is_true())
            
  FreqDisp0 <- c(158, 166)
  FreqDisp <- as.vector(table(Test$TableDisp$ResultDiffDisp))
  expect_that(all.equal(FreqDisp, FreqDisp0), is_true())
  
  FreqRes0 <- c(328,   7,  15)
  FreqRes <- as.vector(table(Test$TableResDisp$ResultDiffResDisp))
  expect_that(all.equal(FreqRes, FreqRes0), is_true())
            
  # Posterior probabilities
            
  ProbMean0 <- c(0.01, 0.01, 0.79, 0.41, 0.35)
  ProbMean <- round(Test$TableMean$ProbDiffMean[1:5], 2)
  expect_that(all.equal(ProbMean, ProbMean0), is_true())
            
  ProbDisp0 <- c(0.81, 0.72, 0.72, 0.72, 0.59)
  ProbDisp <- round(tail(Test$TableDisp$ProbDiffDisp, 5), 2)
  expect_that(all.equal(ProbDisp, ProbDisp0), is_true())
  
  ProbRes0 <- c(0.60, 0.56, 0.63, 1.00, 0.52)
  ProbRes <- round(tail(Test$TableRes$ProbDiffRes, 5), 2)
  expect_that(all.equal(ProbRes, ProbRes0), is_true())
            
  # Log2 fold changes
  Lfc2Mean0 <- c(-0.22, -0.08, -0.92, -0.54, -0.46)
  Lfc2Mean <- round(Test$TableMean$MeanLog2FC[1:5], 2)
  expect_that(all.equal(Lfc2Mean, Lfc2Mean0), is_true())
            
  Lfc2Disp0 <- c(1.19, 0.92, 1.02, 0.86, 0.26)
  Lfc2Disp <- round(tail(Test$TableDisp$DispLog2FC, 5),2)
  expect_that(all.equal(Lfc2Disp, Lfc2Disp0), is_true())
            
})




