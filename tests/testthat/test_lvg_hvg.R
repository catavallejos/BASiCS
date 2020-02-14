context("Basic example of HVG/LVG detection")

test_that("HVG/LVG detection is correct", {
  data(ChainSC)

  DetectHVG <- BASiCS_DetectHVG(ChainSC, PercentileThreshold = NULL,
                                VarThreshold = 0.60,
                                EFDR = 0.10, Plot = FALSE)
  DetectLVG <- BASiCS_DetectLVG(ChainSC, PercentileThreshold = NULL,
                                VarThreshold = 0.40,
                                EFDR = 0.10, Plot = FALSE)
  
  FreqHVG0 <- c(347, 3)
  FreqHVG <- as.vector(table(DetectHVG$Table$HVG))  
  expect_equal(FreqHVG, FreqHVG0)
  
  FreqLVG0 <- c( 119, 231)
  FreqLVG <- as.vector(table(DetectLVG$Table$LVG))  
  expect_equal(FreqLVG, FreqLVG0)

  ProbHVG0 <- c(0.97, 0.96, 0.69, 0.68, 0.65)
  ProbHVG <- round(DetectHVG$Table$Prob[1:5], 2)
  expect_equal(ProbHVG, ProbHVG0)
  
  ProbLVG0 <- c(0.92, 0.92, 0.92, 0.91, 0.91)
  ProbLVG <- round(DetectLVG$Table$Prob[141:145], 2)
  expect_equal(ProbLVG, ProbLVG0)  
})

test_that("HVG/LVG detection using epsilons is correct", {
  data(ChainSCReg)
            
  DetectHVG <- BASiCS_DetectHVG(ChainSCReg, EFDR = 0.10, Plot = FALSE)
  DetectLVG <- BASiCS_DetectLVG(ChainSCReg, EFDR = 0.10, Plot = FALSE)
            
  FreqHVG0 <- c(332, 18)
  FreqHVG <- as.vector(table(DetectHVG$Table$HVG))  
  expect_equal(FreqHVG, FreqHVG0)
            
  FreqLVG0 <- c(340, 10)
  FreqLVG <- as.vector(table(DetectLVG$Table$LVG))  
  expect_equal(FreqLVG, FreqLVG0)
            
  ProbHVG0 <- c(1.00, 1.00, 1.00, 1.00, 0.99)
  ProbHVG <- round(DetectHVG$Table$Prob[1:5], 2)
  expect_equal(ProbHVG, ProbHVG0)
            
  ProbLVG0 <- c(1.00, 1.00, 1.00, 1.00, 0.97)
  ProbLVG <- round(DetectLVG$Table$Prob[1:5], 2)
  expect_equal(ProbLVG, ProbLVG0)
})


test_that("HVG/LVG plotting works", {
  data(ChainSCReg)
            
  DetectHVG <- BASiCS_DetectHVG(ChainSCReg, EFDR = 0.10, Plot = TRUE)
  DetectLVG <- BASiCS_DetectLVG(ChainSCReg, EFDR = 0.10, Plot = TRUE)
  expect_is(DetectHVG$Plots[[1]], "gg")
  expect_is(DetectHVG$Plots[[2]], "gg")
  expect_is(DetectLVG$Plots[[1]], "gg")
  expect_is(DetectLVG$Plots[[2]], "gg")
})
