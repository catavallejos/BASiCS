HiddenBASiCS_MCMC_ExtraArgs <- function(Data,
                                        Burn,
                                        GPar,
                                        Regression,
                                        WithSpikes,
                                        PriorDelta = c("log-normal", "gamma"),
                                        PriorDeltaNum = if (PriorDelta == "gamma") 1 else 2,
                                        k = 12,
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
                                          m = rep(0, k),
                                          V = diag(k),
                                          a.sigma2 = 2,
                                          b.sigma2 = 2,
                                          eta = 5
                                        ),
                                        Start = HiddenBASiCS_MCMC_Start(
                                          Data = Data,
                                          eta = PriorParam$eta,
                                          m = PriorParam$m,
                                          V = PriorParam$V,
                                          a.sigma2 = PriorParam$a.sigma2,
                                          b.sigma2 = PriorParam$b.sigma2,
                                          WithSpikes = WithSpikes
                                        ),
                                        Locations = HiddenFindRBFLocations(
                                          Data = Start$mu0, 
                                          k = k
                                        ),
                                        FixLocations = !missing(Locations),
                                        mintol_mu = 1e-3,
                                        mintol_delta = 1e-3,
                                        mintol_nu = 1e-5,
                                        mintol_theta = 1e-4
                                        ) {

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
    if (k <= 3) {
      stop("The number of basis functions needs to be >= 4.")
    }
  }


  # Validity checks
  if (!(
        # PriorParam$mu.mu > 0 & 
        length(PriorParam$mu.mu) == 1 &
        PriorParam$s2.mu > 0 & length(PriorParam$s2.mu) == 1 &
        PriorParam$s2.delta > 0 & length(PriorParam$s2.delta) == 1 &
        PriorParam$a.delta > 0 & length(PriorParam$a.delta) == 1 &
        PriorParam$b.delta > 0 & length(PriorParam$b.delta) == 1 &
        all(PriorParam$p.phi > 0) & length(PriorParam$p.phi) == GPar$n &
        PriorParam$a.s > 0 & length(PriorParam$a.s) == 1 &
        PriorParam$b.s > 0 & length(PriorParam$b.s) == 1 &
        PriorParam$a.theta > 0 & length(PriorParam$a.theta) == 1 &
        PriorParam$b.theta > 0 & length(PriorParam$b.theta) == 1)) {
    stop("Invalid prior hyper-parameter values.")
  }
  if (Regression) {
    if (!(length(PriorParam$m) == k &
         ncol(PriorParam$V) == k & nrow(PriorParam$V) == k &
         PriorParam$a.sigma2 > 0 & length(PriorParam$a.sigma2) == 1 &
         PriorParam$b.sigma2 > 0 & length(PriorParam$b.sigma2) == 1 )) {
      stop("Invalid prior hyper-parameter values.")
    }
  }

  if (!(AR > 0 & AR < 1 & length(AR) == 1)) {
    stop("Invalid AR value. Recommended value: AR = 0.44.")
  }
  if (!(StopAdapt > 0)) {
    stop("Invalid StopAdapt value.")
  }
  if (!(is.logical(StoreChains) & length(StoreChains) == 1)) {
    stop("Invalid StoreChains value.")
  }
  if (!(is.logical(StoreAdapt) & length(StoreAdapt) == 1)) {
    stop("Invalid StoreAdapt value.")
  }
  if (!(file.info(StoreDir)["isdir"])) {
    stop("Invalid StoreDir value.")
  }
  if (!(PriorDelta %in% c("gamma", "log-normal"))) {
    stop("Invalid PriorDelta value.")
  }
  
  if(!(mintol_mu > 0) | !(mintol_delta > 0) | 
     !(mintol_nu > 0) | !(mintol_theta > 0)) {
    stop("Invalid value for mintol parameters (mu, delta, nu or theta)")
  }

  # Definition of parameters that are specific to the no-spikes case
  NoSpikesParam <- HiddenBASiCS_MCMC_NoSpikesParam(
    GPar$BioCounts,
    ConstrainType,
    StochasticRef, 
    GPar$q.bio, 
    Start$mu0, 
    PriorDelta, 
    ConstrainProp)
  ConstrainGene <- NoSpikesParam$ConstrainGene
  NotConstrainGene <- NoSpikesParam$NotConstrainGene
  Constrain <- NoSpikesParam$Constrain

  RefGenes <- NoSpikesParam$RefGenes
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
    WithSpikes = WithSpikes, 
    StochasticRef = StochasticRef,
    ConstrainType = ConstrainType, 
    ConstrainProp = ConstrainProp,
    k = k, 
    Locations = Locations,
    FixLocations = FixLocations,
    variance = variance,
    Start = Start,
    mintol_mu = mintol_mu, 
    mintol_delta = mintol_delta,
    mintol_nu = mintol_nu, 
    mintol_theta = mintol_theta,
    ConstrainGene = NoSpikesParam$ConstrainGene,
    NotConstrainGene = NoSpikesParam$NotConstrainGene,
    Constrain = NoSpikesParam$Constrain,
    RefGenes = NoSpikesParam$RefGenes,
    RefGene = NoSpikesParam$RefGene,
    Index = Index
  )

}
