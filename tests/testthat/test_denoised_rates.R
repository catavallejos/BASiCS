context("Denoised rates\n")

test_that("Denoised rates match the given seed", 
{
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  set.seed(18)
  Chain <- BASiCS_MCMC(Data, N = 10000, Thin = 10, Burn = 5000, 
                       PrintProgress = FALSE)
   
  DR <- BASiCS_DenoisedRates(Data, Chain)

  # Checks for 2 arbitrary sets of genes / cells
  
  DRcheck0 <- c(2.809, 1.616, 2.499, 4.925, 8.812)
  DRcheck <- as.vector(round(DR[1:5,1], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
            
  DRcheck0 <- c(2.104, 2.922, 3.592, 2.528, 3.282)
  DRcheck <- as.vector(round(DR[10,1:5], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
})
