context("Differential test")

test_that("Differential test is correct", 
          {
            data(ChainSC)
            data(ChainRNA)
            
            Test <- BASiCS_TestDE(Chain1 = ChainSC, Chain2 = ChainRNA,
                                  GroupLabel1 = 'SC', GroupLabel2 = 'P&S',
                                  EpsilonM = log2(1.5), EpsilonD = log2(1.5), 
                                  OffSet = TRUE, Plot = FALSE, PlotOffset = FALSE)
            
            # Classification frequency
            
            FreqMean0 <- c(469,  19,  12)
            FreqMean <- as.vector(table(Test$TableMean$ResultDiffMean))
            expect_that(all.equal(FreqMean, FreqMean0), is_true())
            
            FreqDisp0 <- c(359, 110)
            FreqDisp <- as.vector(table(Test$TableDisp$ResultDiffDisp))
            expect_that(all.equal(FreqDisp, FreqDisp0), is_true())
            
            # Posterior probabilities
            
            ProbMean0 <- c(1.00, 1.00, 1.00, 1.00, 1.00, 
                           1.00, 1.00, 0.99, 0.98 ,0.97)
            ProbMean <- Test$TableMean$ProbDiffMean[1:10]
            expect_that(all.equal(ProbMean, ProbMean0), is_true())
            
            ProbDisp0 <- c(0.38, 0.38, 0.37, 0.37, 0.37,
                           0.35, 0.34, 0.31, 0.29, 0.29)
            ProbDisp <- tail(Test$TableDisp$ProbDiffDisp, 10)
            expect_that(all.equal(ProbDisp, ProbDisp0), is_true())
            
            # Log2 fold changes
            
            Lfc2Mean0 <- c( 4.966, -1.607, -0.802,  1.490, -2.510, 
                           -1.605,  0.799, -0.825, -0.900, -0.958)
            Lfc2Mean <- Test$TableMean$MeanLog2FC[1:10]
            expect_that(all.equal(Lfc2Mean, Lfc2Mean0), is_true())
            
            Lfc2Disp0 <- c(0.028,  0.244,  0.120,  0.285, -0.081, 
                           0.116, -0.104,  0.158, -0.138,  0.287)
            Lfc2Disp <- tail(Test$TableDisp$DispLog2FC, 10)
            expect_that(all.equal(Lfc2Disp, Lfc2Disp0), is_true())
            
          })




