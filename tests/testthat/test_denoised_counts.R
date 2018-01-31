context("Denoised counts\n")

test_that("Denoised counts match the given seed", 
          {
            Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
            set.seed(18)
            Chain <- BASiCS_MCMC(Data, N = 2000, Thin = 10, Burn = 1000, 
                                 PrintProgress = FALSE)
            
            DC <- BASiCS_DenoisedCounts(Data, Chain)
            
            # Checks for 2 arbitrary sets of genes / cells
            
            DCcheck0 <- c(0.000, 0.000, 0.000, 4.672, 4.672)
            DCcheck <- as.vector(round(DC[1:5,1], 3))
            expect_that(all.equal(DCcheck, DCcheck0), is_true())
            
            DCcheck0 <- c(0.000, 1.863, 0.000, 0.000, 2.629)
            DCcheck <- as.vector(round(DC[10,1:5], 3))
            expect_that(all.equal(DCcheck, DCcheck0), is_true())
          })
