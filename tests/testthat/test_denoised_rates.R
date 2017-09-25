context("Denoised rates\n")

test_that("Denoised rates match the given seed", 
{
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE, Example = 1)
  set.seed(18)
  Chain <- BASiCS_MCMC(Data, N = 10000, Thin = 10, Burn = 5000, 
                       PrintProgress = FALSE)
   
  DR <- BASiCS_DenoisedRates(Data, Chain)

  # Checks for 2 arbitrary sets of genes / cells
  
  DRcheck0 <- c(2.711, 1.521, 2.676, 5.123, 9.303)
  DRcheck <- as.vector(round(DR[1:5,1], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
            
  DRcheck0 <- c(2.156, 2.974, 3.547, 2.566, 3.314)
  DRcheck <- as.vector(round(DR[10,1:5], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
})
