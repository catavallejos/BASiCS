context("Output of BASiCS_MCMC\n")

test_that("Valid BASiCS_MCMC output object", 
{
  # Data example: spikes
  set.seed(7)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  # Running the samples: spikes + no regression
  set.seed(18)
  Chain <- run_MCMC(Data, N = 200, Thin = 2, Burn = 100, 
                    Regression = FALSE, WithSpikes = TRUE, 
                    PrintProgress = FALSE, StoreAdapt = TRUE)
  # Checking parameter names
  ParamNames <- c("mu", "delta", "phi", "s", "nu", "theta")
  expect_equal(names(Chain@parameters), ParamNames)
  # Calculating a posterior summary
  PostSummary <- Summary(Chain)  
  expect_equal(names(PostSummary@parameters), ParamNames)
  
  # Running the samples: spikes + regression
  set.seed(18)
  Chain <- run_MCMC(Data, N = 200, Thin = 2, Burn = 100, 
                    Regression = TRUE, WithSpikes = TRUE, 
                    PrintProgress = FALSE, StoreAdapt = TRUE)
  # Checking parameter names
  ParamNames <- c("mu", "delta", "phi", "s", "nu", "theta",
                  "beta", "sigma2", "epsilon", "designMatrix")
  expect_equal(names(Chain@parameters), ParamNames)
  # Calculating a posterior summary
  PostSummary <- Summary(Chain)  
  expect_equal(names(PostSummary@parameters), ParamNames[-length(ParamNames)])

  # Data example: no-spikes
  set.seed(7)
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE)
  # Running the samples: no-spikes + no regression
  set.seed(18)
  Chain <- run_MCMC(Data, N = 200, Thin = 2, Burn = 100, 
                    Regression = FALSE, WithSpikes = FALSE, 
                    PrintProgress = FALSE, StoreAdapt = TRUE)
  # Checking parameter names
  ParamNames <- c("mu", "delta", "s", "nu", "theta", "RefFreq")
  expect_equal(names(Chain@parameters), ParamNames)
  # Calculating a posterior summary
  PostSummary <- Summary(Chain)  
  expect_equal(names(PostSummary@parameters), ParamNames[-length(ParamNames)])
 
  # Running the samples: no-spikes + regression
  set.seed(18)
  Chain <- run_MCMC(Data, N = 200, Thin = 2, Burn = 100, 
                    Regression = TRUE, WithSpikes = FALSE, 
                    PrintProgress = FALSE, StoreAdapt = TRUE)
  # Checking parameter names
  ParamNames <- c("mu", "delta", "s", "nu", "theta",
                  "beta", "sigma2", "epsilon", "RefFreq", "designMatrix")
  expect_equal(names(Chain@parameters), ParamNames)
  # Calculating a posterior summary
  PostSummary <- Summary(Chain)
  ind_names <- setdiff(
    seq_along(ParamNames), 
    c(length(ParamNames)-1, length(ParamNames))
  )
  expect_equal(
    names(PostSummary@parameters), 
    ParamNames[ind_names]
  )
 
})


