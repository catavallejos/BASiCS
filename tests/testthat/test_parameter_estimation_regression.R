context("Parameter estimation (regression case)\n")

test_that("paramater estimations match the given seed", 
          {
            Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
            set.seed(12)
            Chain <- BASiCS_MCMC(Data, N = 2000, Thin = 10, Burn = 1000, 
                                 PrintProgress = FALSE, Regression = TRUE)
            PostSummary <- Summary(Chain)
            
            # Check if parameter estimates match for the first 5 genes and cells
            Mu <- c(7.309,  4.144,  3.783,  5.009, 18.222)
            MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
            expect_that(all.equal(MuObs, Mu), is_true())
            
            Delta <- c(1.126, 2.326, 1.171, 1.565, 0.556)
            DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "delta")[1:5,1],3))
            expect_that(all.equal(DeltaObs, Delta), is_true())
            
            Phi <- c(1.158, 1.008, 0.974, 1.082, 0.954)
            PhiObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "phi")[1:5,1],3))
            expect_that(all.equal(PhiObs, Phi), is_true())
            
            S <- c( 0.243, 0.466, 0.056, 0.195, 0.535)
            SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
            expect_that(all.equal(SObs, S), is_true())
            
            Theta <- 0.291
            ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
            expect_that(all.equal(ThetaObs, Theta), is_true())
          })

