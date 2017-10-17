context("Denoised counts\n")

test_that("Denoised counts match the given seed", 
          {
            Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
            set.seed(18)
            Chain <- BASiCS_MCMC(Data, N = 10000, Thin = 10, Burn = 5000, 
                                 PrintProgress = FALSE)
            
            DC <- BASiCS_DenoisedCounts(Data, Chain)
            
            # Checks for 2 arbitrary sets of genes / cells
            
            DCcheck0 <- c(0.000, 0.000, 0.000, 4.765, 4.765)
            DCcheck <- as.vector(round(DC[1:5,1], 3))
            expect_that(all.equal(DCcheck, DCcheck0), is_true())
            
            DCcheck0 <- c(0.000, 2.200, 0.000, 0.000, 2.638)
            DCcheck <- as.vector(round(DC[10,1:5], 3))
            expect_that(all.equal(DCcheck, DCcheck0), is_true())
          })
