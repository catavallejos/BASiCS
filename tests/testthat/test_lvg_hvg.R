context("Basic example of HVG/LVG detection\n")

test_that("HVG/LVG detection is correct", 
{
  data(ChainSC)

  DetectHVG <- BASiCS_DetectHVG(ChainSC, VarThreshold = 0.60, 
                                EFDR = 0.10, Plot = FALSE)
  DetectLVG <- BASiCS_DetectLVG(ChainSC, VarThreshold = 0.40, 
                                EFDR = 0.10, Plot = FALSE)
  
  FreqHVG0 <- c(347,   3)
  FreqHVG <- as.vector(table(DetectHVG$Table$HVG))  
  expect_that(all.equal(FreqHVG, FreqHVG0), is_true())
  
  FreqLVG0 <- c( 119, 231)
  FreqLVG <- as.vector(table(DetectLVG$Table$LVG))  
  expect_that(all.equal(FreqLVG, FreqLVG0), is_true())         

  ProbHVG0 <- c(0.97, 0.96, 0.69, 0.68, 0.65)
  ProbHVG <- round(DetectHVG$Table$Prob[1:5], 2)
  expect_that(all.equal(ProbHVG, ProbHVG0), is_true())
  
  ProbLVG0 <- c(0.92, 0.92, 0.92, 0.91, 0.91)
  ProbLVG <- round(DetectLVG$Table$Prob[141:145], 2)
  expect_that(all.equal(ProbLVG, ProbLVG0), is_true())
  
})
