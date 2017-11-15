context("Parameter estimation (regression case)\n")

test_that("paramater estimations match the given seed", 
          {
            Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
            set.seed(12)
            Chain <- BASiCS_MCMC(Data, N = 10000, Thin = 10, Burn = 5000, 
                                 PrintProgress = FALSE, Regression = TRUE,
                                 k = 12, Var = 1.2)
            PostSummary <- Summary(Chain)
            
            # Check if parameter estimates match for the first 5 genes and cells
            Mu <- c(7.446,  4.413,  4.009,  4.917, 18.996)
            MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
            expect_that(all.equal(MuObs, Mu), is_true())
            
            Delta <- c(1.163, 2.436, 0.950, 1.609, 0.538)
            DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "delta")[1:5,1],3))
            expect_that(all.equal(DeltaObs, Delta), is_true())
            
            Phi <- c(1.109, 1.018, 0.907, 0.965, 0.913)
            PhiObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "phi")[1:5,1],3))
            expect_that(all.equal(PhiObs, Phi), is_true())
            
            S <- c(0.320, 0.602, 0.086, 0.220, 0.536)
            SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
            expect_that(all.equal(SObs, S), is_true())
            
            Theta <- 0.438
            ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
            expect_that(all.equal(ThetaObs, Theta), is_true())
          })

