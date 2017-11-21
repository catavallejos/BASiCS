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
            Mu <- c(7.852,  4.456,  4.179,  5.009, 18.662)
            MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
            expect_that(all.equal(MuObs, Mu), is_true())
            
            Delta <- c(0.979, 2.281, 1.076, 1.457, 0.518)
            DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "delta")[1:5,1],3))
            expect_that(all.equal(DeltaObs, Delta), is_true())
            
            Phi <- c(1.122, 1.038, 0.872, 1.005, 0.904)
            PhiObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "phi")[1:5,1],3))
            expect_that(all.equal(PhiObs, Phi), is_true())
            
            S <- c(0.290, 0.576, 0.087, 0.224, 0.521)
            SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
            expect_that(all.equal(SObs, S), is_true())
            
            Theta <- 0.456
            ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
            expect_that(all.equal(ThetaObs, Theta), is_true())
          })

