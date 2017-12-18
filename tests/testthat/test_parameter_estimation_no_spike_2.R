context("Parameter estimation (no spike case), original data has spikes\n")

test_that("paramater estimations match the given seed", 
{
  Data1 <- makeExampleBASiCS_Data(WithSpikes = TRUE, 
                                  WithBatch = TRUE)
  Data2 <- newBASiCS_Data(Counts = assay(Data1)[!isSpike(Data1),],
                          Tech = isSpike(Data1)[!isSpike(Data1)],
                          BatchInfo = metadata(Data1)$BatchInfo)
  expect_that(all.equal(assay(Data1)[!isSpike(Data1),], assay(Data2)), is_true())

  set.seed(16)
  Chain1 <- BASiCS_MCMC(Data1, N = 2000, Thin = 10, Burn = 1000, 
                        PrintProgress = FALSE, WithSpikes = FALSE)
  PostSummary1 <- Summary(Chain1)
  set.seed(16)
  Chain2 <- BASiCS_MCMC(Data2, N = 2000, Thin = 10, Burn = 1000, 
                        PrintProgress = FALSE, WithSpikes = FALSE)
  PostSummary2 <- Summary(Chain2)
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu1 <- as.vector(round(displaySummaryBASiCS(PostSummary1, "mu")[,1],3))
  Mu2 <- as.vector(round(displaySummaryBASiCS(PostSummary2, "mu")[,1],3))
  expect_that(all.equal(Mu1, Mu2), is_true())
            
  Delta1 <- as.vector(round(displaySummaryBASiCS(PostSummary1, "delta")[,1],3))
  Delta2 <- as.vector(round(displaySummaryBASiCS(PostSummary2, "delta")[,1],3))
  expect_that(all.equal(Delta1, Delta2), is_true())

  Phi1 <- as.vector(round(displaySummaryBASiCS(PostSummary1, "phi")[,1],3))
  Phi2 <- as.vector(round(displaySummaryBASiCS(PostSummary2, "phi")[,1],3))
  expect_that(all.equal(Phi1, Phi2), is_true())
            
  S1 <- as.vector(round(displaySummaryBASiCS(PostSummary1, "s")[,1],3))
  S2 <- as.vector(round(displaySummaryBASiCS(PostSummary2, "s")[,1],3))
  expect_that(all.equal(S1, S2), is_true())
  
  Theta1 <- as.vector(round(displaySummaryBASiCS(PostSummary1, "theta")[,1],3))
  Theta2 <- as.vector(round(displaySummaryBASiCS(PostSummary2, "theta")[,1],3))
  expect_that(all.equal(Theta1, Theta2), is_true())
})

