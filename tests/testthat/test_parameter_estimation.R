context("Basic example of parameter estimation")

test_that("paramater estimations are correct", {
    Data = makeExampleBASiCS_Data(WithSpikes = TRUE, Case1 = TRUE)
    MCMC_Output <- BASiCS_MCMC(Data, N = 10000, Thin = 10, Burn = 5000, PrintProgress = FALSE)
    MCMC_Summary <- Summary(MCMC_Output, prob = 0.995)
    # Parameters used for the simulation
    Mu =  c( 8.36,  10.65,   4.88,   6.29,  21.72,  12.93,  30.19,  83.92,   3.89,   6.34,
             57.87,  12.45,   8.08,   7.31,  15.56,  15.91,  12.24,  15.96,  19.14,   4.20,
             6.75,  27.74,   8.88,  21.21,  19.89,   7.14,  11.09,   7.19,  20.64,  73.90,
             9.05,   6.13,  16.40,   6.08,  17.89,   6.98,  10.25,  14.05,   8.14,   5.67,
             6.95,  11.16,  11.83,   7.56, 159.05,  16.41,   4.58,  15.46,  10.96,  25.05,
             18.54,   8.50,   4.05,   5.37,   4.63,   4.08,   3.75,   5.02,  27.74,  10.28)
    Delta = c(1.29, 0.88, 1.51, 1.49, 0.54, 0.40, 0.85, 0.27, 0.53, 1.31,
              0.26, 0.81, 0.72, 0.70, 0.96, 0.58, 1.15, 0.82, 0.25, 5.32,
              1.13, 0.31, 0.66, 0.27, 0.76, 1.39, 1.18, 1.57, 0.55, 0.17,
              1.40, 1.47, 0.57, 2.55, 0.62, 0.77, 1.47, 0.91, 1.53, 2.89,
              1.43, 0.77, 1.37, 0.57, 0.15, 0.33, 3.99, 0.47, 0.87, 0.86)
    Phi = c(1.09, 1.16, 1.19, 1.14, 0.87, 1.10, 0.48, 1.06, 0.94, 0.97,
            1.09, 1.16, 1.19, 1.14, 0.87, 1.10, 0.48, 1.06, 0.94, 0.97)
    S = c(0.38, 0.40, 0.38, 0.39, 0.34, 0.39, 0.31, 0.39, 0.40, 0.37,
          0.38, 0.40, 0.38, 0.39, 0.34, 0.39, 0.31, 0.39, 0.40, 0.37)

    s = sum(Mu >= MCMC_Summary@mu[,2] & Mu <= MCMC_Summary@mu[,3]) / 100
    expect_that(s >= 0.9, is_true())

    s = sum(Delta >= MCMC_Summary@delta[,2] & Delta <= MCMC_Summary@delta[,3])
    s = s / 100
    expect_that(s >= 0.9, is_true())

    s = sum(Phi >= MCMC_Summary@phi[,2] & Phi <= MCMC_Summary@phi[,3]) / 20
    expect_that(s >= 0.9, is_true())

    s = sum(S >= MCMC_Summary@s[,2] & S <= MCMC_Summary@s[,3]) / 20  >= 0.9
    expect_that(s >= 0.9, is_true())
})
