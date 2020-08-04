context("Basic example of HVG/LVG detection")

test_that("HVG/LVG consistency of input values= no regression", {
  
  data(ChainSC)
  # If none of the required args is provided. 
  expect_error(BASiCS_DetectVG(ChainSC),
               regexp = "argument \"N\" is missing, with no default")
  
})

test_that("HVG/LVG detection is correct", {
  data(ChainSC)

  DetectHVG <- BASiCS_DetectHVG(ChainSC, PercentileThreshold = NULL,
                                VarThreshold = 0.60,
                                EFDR = 0.10, Plot = FALSE)
  DetectLVG <- BASiCS_DetectLVG(ChainSC, PercentileThreshold = NULL,
                                VarThreshold = 0.40,
                                EFDR = 0.10, Plot = FALSE)
  
  FreqHVG0 <- c(347, 3)
  FreqHVG <- as.vector(table(DetectHVG@Table$HVG))  
  expect_equal(FreqHVG, FreqHVG0)
  
  FreqLVG0 <- c( 119, 231)
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

  DetectHVG <- BASiCS_DetectHVG(ChainSCReg, PercentileThreshold = 0.90,
                                EFDR = 0.10, Plot = FALSE)
  DetectLVG <- BASiCS_DetectLVG(ChainSCReg, PercentileThreshold = 0.10,
                                EFDR = 0.10, Plot = FALSE)

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
  DetectHVG <- BASiCS_DetectHVG(ChainSCReg, EFDR = 0.10, Plot = FALSE)
  expect_is(format(DetectHVG), "data.frame")
  expect_error(DetectHVG[1:10, ], NA)
  p <- 0.8
  d <- format(DetectHVG, ProbThreshold = p)
  expect_true(all(d$Prob > p))
  rd <- rowData(DetectHVG)
  ind <- seq_len(nrow(rd))
  rd$ind <- ind
  rowData(DetectHVG) <- rd
  expect_equal(rowData(DetectHVG)$ind, ind)
})

test_that("HVG/LVG plotting works", {
  data(ChainSCReg)
            
  DetectHVG <- BASiCS_DetectHVG(ChainSCReg, EFDR = 0.10, Plot = TRUE)
  DetectLVG <- BASiCS_DetectLVG(ChainSCReg, EFDR = 0.10, Plot = TRUE)
  expect_is(DetectHVG@Extras$Plots[[1]], "gg")
  expect_is(DetectHVG@Extras$Plots[[2]], "gg")
  expect_is(DetectLVG@Extras$Plots[[1]], "gg")
  expect_is(DetectLVG@Extras$Plots[[2]], "gg")
})



test_that("VarThresholdSearch works", {
  data(ChainSCReg)
  grid <- seq(0.25, 0.75, length.out = 3)
  DetectHVG <- BASiCS_VarThresholdSearchHVG(ChainSCReg,
    VarThresholdsGrid = grid)
  expect_is(DetectHVG, "data.frame")
  DetectLVG <- BASiCS_VarThresholdSearchLVG(ChainSCReg,
    VarThresholdsGrid = grid)
  expect_is(DetectHVG, "data.frame")
})
