context("Basic example of HVG/LVG detection")

test_that("HVG/LVG - check if key parameters are missing", {
  
  data(ChainSC)
  expect_error(BASiCS_DetectVG(ChainSC, Task = "HVG"),
               regexp = "A value must be provided for")
  expect_error(BASiCS_DetectVG(ChainSC, Task = "HVG",
                               PercentileThreshold = 0.9),
               regexp = "\'Chain\' does not include residual")
  expect_error(BASiCS_DetectVG(ChainSC, Task = "HVG",
                               VarThreshold = 0.6), NA)
  expect_message(BASiCS_DetectVG(ChainSC, Task = "HVG",
                                 VarThreshold = 0.6,
                                 EFDR = NULL), 
               regexp = "EFDR = NULL for")

  data(ChainSCReg)
  expect_error(BASiCS_DetectVG(ChainSCReg, Task = "HVG"),
               regexp = "A value must be provided for")
  expect_error(BASiCS_DetectVG(ChainSCReg, Task = "HVG",
                               VarThreshold = 0.9),
               regexp = "\'Chain\' includes residual")
  expect_error(BASiCS_DetectVG(ChainSCReg, Task = "HVG",
                               PercentileThreshold = 0.9), NA)
  expect_message(BASiCS_DetectVG(ChainSCReg, Task = "HVG",
                                 PercentileThreshold = 0.9,
                                 EFDR = NULL), 
                 regexp = "EFDR = NULL for")
  
})

test_that("HVG/LVG detection is correct", {
  data(ChainSC)

  DetectHVG <- BASiCS_DetectHVG(ChainSC, VarThreshold = 0.60)
  DetectLVG <- BASiCS_DetectLVG(ChainSC, VarThreshold = 0.40)
  
  FreqHVG0 <- c(347, 3)
  FreqHVG <- as.vector(table(DetectHVG@Table$HVG))  
  expect_equal(FreqHVG, FreqHVG0)
  
  FreqLVG0 <- c(146, 204)
  FreqLVG <- as.vector(table(DetectLVG@Table$LVG))  
  expect_equal(FreqLVG, FreqLVG0)

  ProbHVG0 <- c(0.97, 0.96, 0.69, 0.68, 0.65)
  ProbHVG <- round(DetectHVG@Table$Prob[1:5], 2)
  expect_equal(ProbHVG, ProbHVG0)
  
  ProbLVG0 <- c(0.92, 0.92, 0.92, 0.91, 0.91)
  ProbLVG <- round(DetectLVG@Table$Prob[141:145], 2)
  expect_equal(ProbLVG, ProbLVG0)

  expect_true(!is.null(DetectHVG@Table$Delta))
  expect_true(!is.null(DetectHVG@Table$Sigma))
})

test_that("HVG/LVG detection using epsilons is correct", {
  data(ChainSCReg)

  DetectHVG <- BASiCS_DetectHVG(ChainSCReg, PercentileThreshold = 0.90)
  DetectLVG <- BASiCS_DetectLVG(ChainSCReg, PercentileThreshold = 0.10)

  FreqHVG0 <- c(332, 18)
  FreqHVG <- as.vector(table(DetectHVG@Table$HVG))  
  expect_equal(FreqHVG, FreqHVG0)
            
  FreqLVG0 <- c(340, 10)
  FreqLVG <- as.vector(table(DetectLVG@Table$LVG))  
  expect_equal(FreqLVG, FreqLVG0)
            
  ProbHVG0 <- c(1.00, 1.00, 1.00, 1.00, 0.99)
  ProbHVG <- round(DetectHVG@Table$Prob[1:5], 2)
  expect_equal(ProbHVG, ProbHVG0)
            
  ProbLVG0 <- c(1.00, 1.00, 1.00, 1.00, 0.97)
  ProbLVG <- round(DetectLVG@Table$Prob[1:5], 2)
  expect_equal(ProbLVG, ProbLVG0)
  
  expect_true(!is.null(DetectHVG@Table$Epsilon))
})

test_that("HVG/LVG utils work", {
  data(ChainSCReg)
  DetectHVG <- BASiCS_DetectHVG(ChainSCReg, PercentileThreshold = 0.9)
  expect_warning(x <- format(DetectHVG))
  expect_is(x, "data.frame")
  expect_error(DetectHVG[1:10, ], NA)
  p <- 0.8
  d <- as.data.frame(DetectHVG, ProbThreshold = p)
  expect_true(all(as.numeric(d$Prob) > p))
  rd <- rowData(DetectHVG)
  ind <- seq_len(nrow(rd))
  rd$ind <- ind
  rowData(DetectHVG) <- rd
  expect_equal(rowData(DetectHVG)$ind, ind)
})

test_that("HVG/LVG plotting works", {
  data(ChainSCReg)
            
  DetectHVG <- BASiCS_DetectHVG(ChainSCReg, PercentileThreshold = 0.9,
                                Plot = TRUE)
  DetectLVG <- BASiCS_DetectLVG(ChainSCReg, PercentileThreshold = 0.9,
                                Plot = TRUE)
  expect_is(DetectHVG@Extras$Plots[[1]], "gg")
  expect_is(DetectHVG@Extras$Plots[[2]], "gg")
  expect_is(DetectLVG@Extras$Plots[[1]], "gg")
  expect_is(DetectLVG@Extras$Plots[[2]], "gg")
})



test_that("VarThresholdSearch works", {
  data(ChainSC)
  grid <- seq(0.25, 0.75, length.out = 3)
  DetectHVG <- BASiCS_VarThresholdSearchHVG(ChainSC,
    VarThresholdsGrid = grid)
  expect_is(DetectHVG, "data.frame")
  DetectLVG <- BASiCS_VarThresholdSearchLVG(ChainSC,
    VarThresholdsGrid = grid)
  expect_is(DetectHVG, "data.frame")
})
