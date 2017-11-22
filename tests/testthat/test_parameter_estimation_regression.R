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
            Mu <- c(7.826,  4.230,  4.174,  4.994, 19.192)
            MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
            expect_that(all.equal(MuObs, Mu), is_true())
            
            Delta <- c(0.946, 2.558, 0.842, 1.469, 0.551)
            DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "delta")[1:5,1],3))
            expect_that(all.equal(DeltaObs, Delta), is_true())
            
            Phi <- c(1.115, 1.032, 0.925, 0.979, 0.922)
            PhiObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "phi")[1:5,1],3))
            expect_that(all.equal(PhiObs, Phi), is_true())
            
            S <- c(0.280, 0.584, 0.078, 0.205, 0.511)
            SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
            expect_that(all.equal(SObs, S), is_true())
            
            Theta <- 0.417
            ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
            expect_that(all.equal(ThetaObs, Theta), is_true())
          })

