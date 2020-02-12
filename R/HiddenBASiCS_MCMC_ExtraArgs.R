HiddenBASiCS_MCMC_ExtraArgs <- function(
    Data,
    Burn,
    GPar,
    Regression,
    WithSpikes,
    PriorDelta = c("log-normal", "gamma"),
    PriorDeltaNum = if (PriorDelta == "gamma") 1 else 2,
    k = 10,
    ## Duplicated arg here for backwards-compatibility.
    ## If both specified, Var is ignored.
    Var = 1.2,
    variance = Var,
    StochasticRef = TRUE,
    ConstrainType = 1,
    ConstrainProp = 0.2,
    AR = 0.44,
    StopAdapt = Burn,
    StoreChains = FALSE,
    StoreAdapt = FALSE,
    StoreDir = getwd(),
    RunName = "",
    PrintProgress = TRUE,
    MinGenesPerRBF = 100,
    PriorParam = list(
      mu.mu = 0,
      s2.mu = 0.5,
      s2.delta = 0.5,
      a.delta = 1,
      b.delta = 1,
      p.phi = rep(1, times = GPar$n),
      a.s = 1,
      b.s = 1,
      a.theta = 1,
      b.theta = 1,
      RBFNTiles = TRUE,
      FixLocations = FALSE,
      m = numeric(k),
      V = diag(k),
      a.sigma2 = 1,
      b.sigma2 = 1,
      eta = eta,
      locations = numeric(k - 2)
    ),
    eta = 5,
    Start = HiddenBASiCS_MCMC_Start(
      Data = Data,
      PriorParam = PriorParam,
      WithSpikes = WithSpikes,
      Regression = Regression
    ),
    GeneExponent = 1,
    CellExponent = 1,
    mintol_mu = 1e-3,
    mintol_delta = 1e-3,
    mintol_nu = 1e-5,
    mintol_theta = 1e-4
  ) {

  .stop_k(k)
  lm <- log(Start$mu0)
  locations <- .estimateRBFLocations(lm, k)
  if (!is.na(MinGenesPerRBF)) {
    d <- (locations[[2]] - locations[[1]]) / 2
    retain <- vapply(
      locations,
      function(location) {
        sum(lm > location - d & lm < location + d) > MinGenesPerRBF
      },
      logical(1)
    )
    # browser()
    # plot(log(Start$mu0), log(Start$delta0))
    # plot_distn <- function(mean, sd) {
    #   x <- seq(0, 10, length=200)
    #   y <- dnorm(x, mean=mean, sd=sd)
    #   lines(x, y, type="l", lwd=2)
    # }
    # abline(v = locations)
    # abline(v = locations + d, col = "grey80", lty = "dashed")
    # abline(v = locations - d, col = "grey80", lty = "dashed")
    # tmp <- sapply(locations, function(location) plot_distn(mean = location, sd = 0.5))
    locations <- locations[retain]
  }
  k <- length(locations) + 2

  PriorDelta <- match.arg(PriorDelta)
  if (missing(PriorDelta) & !Regression) {
    message(
      "-------------------------------------------------------------\n",
      "NOTE: default choice PriorDelta = 'log-normal'  (recommended value). \n",
      "Vallejos et al (2015) used a 'gamma' prior instead.\n",
      "-------------------------------------------------------------\n"
    )
  }
  if (Regression) {
    .stop_k(k)
  }

  if (Regression) {
    PriorParam$m <- rep(0, k); PriorParam$V <- diag(k)
    PriorParam$a.sigma2 <- 2; PriorParam$b.sigma2 <- 2
    PriorParam$eta <- eta
  }

  # Validity checks
  assertthat::assert_that(
    length(PriorParam$mu.mu) == 1,
    PriorParam$s2.mu > 0,
    length(PriorParam$s2.mu) == 1,
    PriorParam$s2.delta > 0,
    length(PriorParam$s2.delta) == 1,
    PriorParam$a.delta > 0,
    length(PriorParam$a.delta) == 1,
    PriorParam$b.delta > 0,
    length(PriorParam$b.delta) == 1,
    all(PriorParam$p.phi > 0),
    length(PriorParam$p.phi) == GPar$n,
    PriorParam$a.s > 0,
    length(PriorParam$a.s) == 1,
    PriorParam$b.s > 0,
    length(PriorParam$b.s) == 1,
    PriorParam$a.theta > 0,
    length(PriorParam$a.theta) == 1,
    PriorParam$b.theta > 0,
    length(PriorParam$b.theta) == 1,
    GeneExponent > 0,
    length(GeneExponent) == 1,
    CellExponent > 0,
    length(CellExponent) == 1,
    StopAdapt > 0,
    length(StoreChains) == 1,
    length(StoreAdapt) == 1,
    is.dir(StoreDir),
    PriorDelta %in% c("gamma", "log-normal"),
    mintol_mu > 0,
    mintol_delta > 0,
    mintol_nu > 0,
    mintol_theta > 0
  )

  if (Regression) {
    assertthat::assert_that(
      length(PriorParam$m) == k,
      ncol(PriorParam$V) == k,
      nrow(PriorParam$V) == k,
      PriorParam$a.sigma2 > 0,
      length(PriorParam$a.sigma2) == 1,
      PriorParam$b.sigma2 > 0,
      length(PriorParam$b.sigma2) == 1
    )
  }

  if (!(AR > 0 & AR < 1 & length(AR) == 1)) {
    stop("Invalid AR value. Recommended value: AR = 0.44.")
  }

  # Definition of parameters that are specific to the no-spikes case
  NoSpikesParam <- HiddenBASiCS_MCMC_NoSpikesParam(
    GPar$BioCounts,
    ConstrainType,
    StochasticRef, 
    GPar$q.bio, 
    Start$mu0, 
    PriorDelta, 
    ConstrainProp
  )
  ConstrainGene <- NoSpikesParam$ConstrainGene
  NotConstrainGene <- NoSpikesParam$NotConstrainGene
  Constrain <- NoSpikesParam$Constrain
  RefGenes <- NoSpikesParam$RefGenes; 
  RefGene <- NoSpikesParam$RefGene
  Index <- seq_len(GPar$q.bio) - 1

  list(
    AR = AR,
    StopAdapt = StopAdapt,
    StoreChains = StoreChains,
    StoreAdapt = StoreAdapt,
    StoreDir = StoreDir,
    RunName = RunName,
    PrintProgress = PrintProgress,
    PriorParam = PriorParam,
    PriorDeltaNum = PriorDeltaNum,
    PriorDelta = PriorDelta,
    StochasticRef = StochasticRef,
    ConstrainType = ConstrainType,
    ConstrainProp = ConstrainProp,
    k = k,
    locations = locations,
    variance = variance,
    Start = Start,
    mintol_mu = mintol_mu,
    mintol_delta = mintol_delta,
    mintol_nu = mintol_nu,
    mintol_theta = mintol_theta,
    ConstrainGene = ConstrainGene,
    NotConstrainGene = NotConstrainGene,
    Constrain = Constrain,
    RefGenes = RefGenes,
    RefGene = RefGene,
    Index = Index,
    GeneExponent = GeneExponent,
    CellExponent = CellExponent
  )
}

## This condition has to be re-used
.stop_k <- function(k) {
  if (k <= 3) {
    stop(
      paste(
        "The number of basis functions needs to be >= 4.",
        "Consider setting MinGenesPerRBF to NA or a lower positive integer.",
        sep = "\n"
      )
    )
  }      
}
