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
            
  FreqMean0 <- c(335,   9,   6)
  FreqMean <- as.vector(table(Test$TableMean$ResultDiffMean))
<<<<<<< HEAD
  expect_true(all.equal(FreqMean, FreqMean0))
=======
  expect_equal(FreqMean, FreqMean0)
>>>>>>> 5d23e114161e35201f17adb88647424ac38c9273
            
  FreqDisp0 <- c(15, 236,  99)
  FreqDisp <- as.vector(table(Test$TableDisp$ResultDiffDisp))
<<<<<<< HEAD
  expect_true(all.equal(FreqDisp, FreqDisp0))
=======
  expect_equal(FreqDisp, FreqDisp0)
>>>>>>> 5d23e114161e35201f17adb88647424ac38c9273
  
  FreqRes0 <- c(338,   4,   8)
  FreqRes <- as.vector(table(Test$TableResDisp$ResultDiffResDisp))
<<<<<<< HEAD
  expect_true(all.equal(FreqRes, FreqRes0))
=======
  expect_equal(FreqRes, FreqRes0)
>>>>>>> 5d23e114161e35201f17adb88647424ac38c9273
            
  # Posterior probabilities
            
  ProbMean0 <- c(0.01, 0.01, 0.79, 0.41, 0.35)
  ProbMean <- round(Test$TableMean$ProbDiffMean[1:5], 2)
<<<<<<< HEAD
  expect_true(all.equal(ProbMean, ProbMean0))
            
  ProbDisp0 <- c(0.72, 0.72, 0.72, 1.00, 0.59)
  ProbDisp <- round(tail(Test$TableDisp$ProbDiffDisp, 5), 2)
  expect_true(all.equal(ProbDisp, ProbDisp0))
  
  ProbRes0 <- c(0.60, 0.56, 0.63, 1.00, 0.52)
  ProbRes <- round(tail(Test$TableRes$ProbDiffRes, 5), 2)
  expect_true(all.equal(ProbRes, ProbRes0))
=======
  expect_equal(ProbMean, ProbMean0)
            
  ProbDisp0 <- c(0.72, 0.72, 0.72, 1.00, 0.59)
  ProbDisp <- round(tail(Test$TableDisp$ProbDiffDisp, 5), 2)
  expect_equal(ProbDisp, ProbDisp0)
  
  ProbRes0 <- c(0.60, 0.56, 0.63, 1.00, 0.52)
  ProbRes <- round(tail(Test$TableRes$ProbDiffRes, 5), 2)
  expect_equal(ProbRes, ProbRes0)
>>>>>>> 5d23e114161e35201f17adb88647424ac38c9273
            
  # Log2 fold changes
  Lfc2Mean0 <- c(-0.22, -0.08, -0.92, -0.54, -0.46)
  Lfc2Mean <- round(Test$TableMean$MeanLog2FC[1:5], 2)
<<<<<<< HEAD
  expect_true(all.equal(Lfc2Mean, Lfc2Mean0))
            
  Lfc2Disp0 <- c(0.92, 1.02, 0.86, 6.23, 0.26)
  Lfc2Disp <- round(tail(Test$TableDisp$DispLog2FC, 5),2)
  expect_true(all.equal(Lfc2Disp, Lfc2Disp0))
=======
  expect_equal(Lfc2Mean, Lfc2Mean0)
            
  Lfc2Disp0 <- c(0.92, 1.02, 0.86, 6.23, 0.26)
  Lfc2Disp <- round(tail(Test$TableDisp$DispLog2FC, 5),2)
  expect_equal(Lfc2Disp, Lfc2Disp0)
>>>>>>> 5d23e114161e35201f17adb88647424ac38c9273
            
})




