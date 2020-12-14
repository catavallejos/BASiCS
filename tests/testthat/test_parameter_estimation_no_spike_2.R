context("Parameter estimation (no-spikes), original has spikes")

test_that("Estimates match (no-spikes)", {
  set.seed(12)
  Data1 <- makeExampleBASiCS_Data(WithSpikes = TRUE, 
                                  WithBatch = TRUE)
  Data2 <- Data1
  altExp(Data2) <- NULL
  rowData(altExp(Data1)) <- NULL
  expect_equal(counts(Data1), counts(Data2))

  set.seed(16)
  Chain1 <- run_MCMC(Data1,
                     N = 1000, Thin = 10, Burn = 500, 
                     Regression = FALSE, WithSpikes = FALSE,
                     PrintProgress = FALSE)
  PostSummary1 <- Summary(Chain1)
  set.seed(16)
  Chain2 <- run_MCMC(Data2, 
                     N = 1000, Thin = 10, Burn = 500,
                     Regression = FALSE, WithSpikes = FALSE, 
                     PrintProgress = FALSE)
  PostSummary2 <- Summary(Chain2)
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu1 <- as.vector(round(displaySummaryBASiCS(PostSummary1, "mu")[,1],3))
  Mu2 <- as.vector(round(displaySummaryBASiCS(PostSummary2, "mu")[,1],3))
  expect_equal(Mu1, Mu2)
            
  Delta1 <- as.vector(round(displaySummaryBASiCS(PostSummary1, "delta")[,1],3))
  Delta2 <- as.vector(round(displaySummaryBASiCS(PostSummary2, "delta")[,1],3))
  expect_equal(Delta1, Delta2)
            
  S1 <- as.vector(round(displaySummaryBASiCS(PostSummary1, "s")[,1],3))
  S2 <- as.vector(round(displaySummaryBASiCS(PostSummary2, "s")[,1],3))
  expect_equal(S1, S2)
  
  Theta1 <- as.vector(round(displaySummaryBASiCS(PostSummary1, "theta")[,1],3))
  Theta2 <- as.vector(round(displaySummaryBASiCS(PostSummary2, "theta")[,1],3))
  expect_equal(Theta1, Theta2)
})

