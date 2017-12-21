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
            
  FreqMean0 <- c(370,  21,   9)
  FreqMean <- as.vector(table(Test$TableMean$ResultDiffMean))
  expect_that(all.equal(FreqMean, FreqMean0), is_true())
            
  FreqDisp0 <- c(180, 190)
  FreqDisp <- as.vector(table(Test$TableDisp$ResultDiffDisp))
  expect_that(all.equal(FreqDisp, FreqDisp0), is_true())
  
  FreqRes0 <- c(374,   9,  17)
  FreqRes <- as.vector(table(Test$TableResDisp$ResultDiffResDisp))
  expect_that(all.equal(FreqRes, FreqRes0), is_true())
            
  # Posterior probabilities
            
  ProbMean0 <- c(1.00, 1.00, 1.00, 1.00, 1.00, 
                 0.99, 0.99, 0.99, 0.97, 0.97)
  ProbMean <- round(Test$TableMean$ProbDiffMean[1:10], 2)
  expect_that(all.equal(ProbMean, ProbMean0), is_true())
            
  ProbDisp0 <- c( 0.51, 0.51, 0.51, 0.51, 0.51, 
                  0.49, 0.49, 0.48, 0.36, 0.35)
  ProbDisp <- round(tail(Test$TableDisp$ProbDiffDisp, 10), 2)
  expect_that(all.equal(ProbDisp, ProbDisp0), is_true())
  
  ProbRes0 <- c( 0.45, 0.45, 0.45, 0.45, 0.43, 
                 0.41, 0.41, 0.40, 0.37, 0.31)
  ProbRes <- round(tail(Test$TableRes$ProbDiffRes, 10), 2)
  expect_that(all.equal(ProbRes, ProbRes0), is_true())
            
  # Log2 fold changes
            
  Lfc2Mean0 <- c(-1.40, -4.28,  1.43, -1.19,  0.86, 
                 -1.21,  0.89, -1.43, -0.87,  0.70)
  Lfc2Mean <- round(Test$TableMean$MeanLog2FC[1:10], 2)
  expect_that(all.equal(Lfc2Mean, Lfc2Mean0), is_true())
            
  Lfc2Disp0 <- c( 0.47,  0.15,  0.21,  0.31, -0.23,  
                  0.24,  0.35, -0.23,  0.24,  0.13)
  Lfc2Disp <- round(tail(Test$TableDisp$DispLog2FC, 10),2)
  expect_that(all.equal(Lfc2Disp, Lfc2Disp0), is_true())
            
})




