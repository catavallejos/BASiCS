context("Denoised rates\n")

test_that("Denoised rates match the given seed", 
{
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  set.seed(18)
  Chain <- BASiCS_MCMC(Data, N = 2000, Thin = 10, Burn = 1000, 
                       PrintProgress = FALSE)
   
  DR <- BASiCS_DenoisedRates(Data, Chain)

  # Checks for 2 arbitrary sets of genes / cells
  
  DRcheck0 <- c(2.769, 1.611, 2.425, 4.684, 8.105)
  DRcheck <- as.vector(round(DR[1:5,1], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
            
  DRcheck0 <- c(2.160, 2.728, 3.847, 2.888, 3.494)
  DRcheck <- as.vector(round(DR[10,1:5], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
})
