context("Storing and loading chains")

test_that("Valid BASiCS_MCMC output object", {
  # Data example: spikes
  set.seed(7)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  # Running the samples: spikes + no regression
  set.seed(18)
  Chain <- run_MCMC(
    Data, N = 20, Thin = 2, Burn = 10, 
    Regression = FALSE, WithSpikes = TRUE, 
    PrintProgress = FALSE, StoreAdapt = TRUE,
    StoreChains = TRUE)
  expect_equal(Chain, readRDS("chain_.Rds"))
  expect_equal(Chain, BASiCS_LoadChain(""))
  # Checking parameter names
  ParamNames <- c("mu", "delta", "phi", "s", "nu", "theta")
  expect_equal(names(Chain@parameters), ParamNames)
  # Calculating a posterior summary
  PostSummary <- Summary(Chain)  
  expect_equal(names(PostSummary@parameters), ParamNames)
  
  # Running the samples: spikes + regression
  set.seed(18)
  Chain <- run_MCMC(
    Data, N = 20, Thin = 2, Burn = 10, 
    Regression = TRUE, WithSpikes = TRUE,
    PrintProgress = FALSE, StoreAdapt = TRUE,
    StoreChains = TRUE)
  expect_equal(Chain, readRDS("chain_.Rds"))
  expect_equal(Chain, BASiCS_LoadChain(""))
  # Checking parameter names
  ParamNames <- c("mu", "delta", "phi", "s", "nu", "theta",
                  "beta", "sigma2", "epsilon", "designMatrix", "RBFLocations")
  expect_equal(names(Chain@parameters), ParamNames)
  # Calculating a posterior summary
  PostSummary <- Summary(Chain)
  expect_equal(names(PostSummary@parameters), ParamNames[-(10:11)])

  # Data example: no-spikes
  set.seed(7)
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE)
  # Running the samples: no-spikes + no regression
  set.seed(18)
  Chain <- run_MCMC(
    Data, N = 20, Thin = 2, Burn = 10, 
    Regression = FALSE, WithSpikes = FALSE, 
    PrintProgress = FALSE, StoreAdapt = TRUE,
    StoreChains = TRUE)
  expect_equal(Chain, readRDS("chain_.Rds"))
  expect_equal(Chain, BASiCS_LoadChain(""))
  # Checking parameter names
  ParamNames <- c("mu", "delta", "s", "nu", "theta", "RefFreq")
  expect_equal(names(Chain@parameters), ParamNames)
  # Calculating a posterior summary
  PostSummary <- Summary(Chain)  
  expect_equal(names(PostSummary@parameters), ParamNames[-length(ParamNames)])
  
  # Running the samples: no-spikes + regression
  set.seed(18)
  Chain <- run_MCMC(
    Data, N = 20, Thin = 2, Burn = 10, 
    Regression = TRUE, WithSpikes = FALSE, 
    PrintProgress = FALSE, StoreAdapt = TRUE,
    StoreChains = TRUE)
  expect_equal(Chain, readRDS("chain_.Rds"))
  expect_equal(Chain, BASiCS_LoadChain(""))
  # Checking parameter names
  ParamNames <- c("mu", "delta", "s", "nu", "theta",
                  "beta", "sigma2", "epsilon", "designMatrix", "RBFLocations", "RefFreq")
  expect_equal(names(Chain@parameters), ParamNames)
  # Calculating a posterior summary
  PostSummary <- Summary(Chain)
  ind_names <- setdiff(
    seq_along(ParamNames), 
    c(length(ParamNames) - 2, length(ParamNames) - 1, length(ParamNames))
  )
  expect_equal(
    names(PostSummary@parameters), 
    ParamNames[ind_names]
  ) 
  file.remove("chain_.Rds")
  file.remove("TableRef_.txt")
})


test_that("reading old files", {
  set.seed(7)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  # Running the samples: spikes + no regression
  set.seed(18)
  Chain <- run_MCMC(
    Data, N = 20, Thin = 2, Burn = 10,
    Regression = FALSE, WithSpikes = TRUE)
  lapply(c("mu", "delta", "nu", "s", "theta", "phi"),
    function(n) {
      write.table(
        Chain@parameters[[n]],
        file = paste0("chain_", n, "_.txt"),
        sep = " "
      )
    }
  )
  expect_equal(Chain, BASiCS_LoadChain())
  lapply(c("mu", "delta", "nu", "s", "theta", "phi"), function(n) {
    file.remove(paste0("chain_", n, "_.txt"))
  })
})

test_that("updating chains works", {
  set.seed(7)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  # Running the samples: spikes + no regression
  set.seed(18)
  Chain <- run_MCMC(
    Data, N = 20, Thin = 2, Burn = 10,
    Regression = FALSE, WithSpikes = TRUE)

  setClass(
    "oldChain",
    contains = NULL,
    representation = representation(
      mu = "matrix",
      delta = "matrix",
      nu = "matrix",
      s = "matrix",
      theta = "matrix",
      phi = "matrix"
    )
  )
  old <- new("oldChain",
    mu = Chain@parameters$mu,
    delta = Chain@parameters$delta,
    nu = Chain@parameters$nu,
    s = Chain@parameters$s,
    theta = Chain@parameters$theta,
    phi = Chain@parameters$phi
  )
  expect_equal(BASiCS:::.updateObject(old), Chain)
  saveRDS(old, "chain_.Rds")
  expect_equal(Chain, BASiCS_LoadChain())
  file.remove("chain_.Rds")
})
