context("Basic example of HVG/LVG detection\n")

test_that("HVG/LVG detection is correct", 
{
  data(ChainSC)

  DetectHVG <- BASiCS_DetectHVG(ChainSC, VarThreshold = 0.60, 
                                EFDR = 0.10, Plot = FALSE)
  DetectLVG <- BASiCS_DetectLVG(ChainSC, VarThreshold = 0.40, 
                                EFDR = 0.10, Plot = FALSE)
  
  FreqHVG0 <- c(493, 7)
  FreqHVG <- as.vector(table(DetectHVG$Table$HVG))  
  expect_that(all.equal(FreqHVG, FreqHVG0), is_true())
  
  FreqLVG0 <- c(224, 276)
  FreqLVG <- as.vector(table(DetectLVG$Table$LVG))  
  expect_that(all.equal(FreqLVG, FreqLVG0), is_true())         

  ProbHVG0 <- c(1.00, 0.96, 0.91, 0.91, 0.89,
                0.88, 0.80, 0.77, 0.77, 0.74)
  ProbHVG <- DetectHVG$Table$Prob[1:10]
  expect_that(all.equal(ProbHVG, ProbHVG0), is_true())
  
  ProbLVG0 <- c(0.98, 0.98, 0.97, 0.97, 0.97,
                0.96, 0.96, 0.96, 0.96, 0.96)
  ProbLVG <- DetectLVG$Table$Prob[141:150]
  expect_that(all.equal(ProbLVG, ProbLVG0), is_true())
  
})
