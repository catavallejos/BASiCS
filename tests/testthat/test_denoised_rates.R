context("Denoised rates\n")

test_that("Denoised rates match the given seed", 
{
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  set.seed(18)
  Chain <- BASiCS_MCMC(Data, N = 10000, Thin = 10, Burn = 5000, 
                       PrintProgress = FALSE)
   
  DR <- BASiCS_DenoisedRates(Data, Chain)

  # Checks for 2 arbitrary sets of genes / cells
  
  DRcheck0 <- c(2.713, 1.550, 2.479, 5.030, 8.787)
  DRcheck <- as.vector(round(DR[1:5,1], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
            
  DRcheck0 <- c(2.160, 2.958, 3.605, 2.586, 3.342)
  DRcheck <- as.vector(round(DR[10,1:5], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
})
