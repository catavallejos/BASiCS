context("Denoised rates\n")

test_that("Denoised rates match the given seed", 
{
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  set.seed(18)
  Chain <- BASiCS_MCMC(Data, N = 2000, Thin = 10, Burn = 1000, 
                       PrintProgress = FALSE)
   
  DR <- BASiCS_DenoisedRates(Data, Chain)

  # Checks for 2 arbitrary sets of genes / cells
  
  DRcheck0 <- c(3.000, 1.752, 2.753, 5.400, 9.464)
  DRcheck <- as.vector(round(DR[1:5,1], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
            
  DRcheck0 <- c(2.088, 2.590, 3.308, 2.358, 3.110)
  DRcheck <- as.vector(round(DR[10,1:5], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
})
