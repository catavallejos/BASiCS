.BASiCS_MCMC_InputCheck <- function(Data, N, Thin,
                                    Burn, Regression, WithSpikes, PriorMu) {
  if (!is(Data, "SingleCellExperiment")) {
    stop("'Data' is not a SingleCellExperiment class object.\n")
  }

  # The following checks are only relevant when the input data was
  # not created using the `newBASiCS_Data` function

  # Check if `Data` contains more than one type of types
  if ((length(SingleCellExperiment::altExpNames(Data)) > 1)) {
    stop(
      "'Data' contains multiple 'altExp' assays but only one is allowed. \n",
      "The latter should include information for spike-in genes. \n",
      "Please remove unwanted the remaining elements. See help(altExp). \n"
    )
  }

  if (WithSpikes & length(altExpNames(Data)) > 0) {
    message(
      "altExp '", altExpNames(Data), "' is assumed to contain spike-in genes.\n",
      "see help(altExp) for details. \n"
    )
  }

  # If SpikeInput slot is missing and WithSpikes == TRUE
  if ((length(altExpNames(Data)) > 0) & WithSpikes &
    is.null(metadata(Data)$SpikeInput)) {
    stop(
      "'Data' does not contain 'SpikeInput' as metadata. \n",
      "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation\n"
    )
  }

  # If isSpike slot is missing and WithSpikes == TRUE
  if ((length(altExpNames(Data)) == 0) & WithSpikes) {
    stop(
      "'Data' does not contain information about spike-in genes \n",
      "Please indicate include this information using 'altExp' \n",
      "or set 'WithSpikes = FALSE' \n.",
      "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation\n"
    )
  }

  # # If isSpike slot does not contain spikes and WithSpikes == TRUE
  # if(sum(Spikes) == 0 & WithSpikes)
  #   stop("'isSpike(Data)' does not contain TRUE values, meaning the sce object \n",
  #        "does not contain spikes. Please indicate in 'isSpike(Data)' which \n",
  #        "genes are spike-ins or set 'WithSpikes = FALSE' \n.",
  #        "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation\n")

  # If BatchInfo slot is missing and WithSpikes == FALSE
  if (!WithSpikes & is.null(colData(Data)$BatchInfo)) {
    stop(
      "'Data' does not contain a BatchInfo vector needed when 'WithSpikes = FALSE'. \n",
      "Please assign the batch information to: 'colData(Data)$BatchInfo = BatchInfo'. \n",
      "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation\n"
    )
  }

  # Checking how counts are stored
  if (!("counts" %in% assayNames(Data))) {
    stop(
      "'Data' does not contain a 'counts' slot. \n",
      "Please make sure to include the raw data in the \n",
      "SingleCellExperiment object under the name 'counts'. \n",
      "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation\n"
    )
  }

  # Further checks on input
  errors <- HiddenChecksBASiCS_Data(Data, WithSpikes)

  if (length(errors) > 0) stop(errors)

  if (!(length(N) == 1 | length(Thin) == 1 | length(Burn) == 1)) {
    stop("Invalid parameter values.\n")
  }
  if (!(N %% Thin == 0 & N >= max(4, Thin))) {
    stop("Please use an integer value for N (N>=4); multiple of thin.\n")
  }
  if (!(Thin %% 1 == 0 & Thin >= 2)) {
    stop("Please use an integer value for Thin (Thin>=2).\n")
  }
  if (!(Burn %% Thin == 0 & Burn < N & Burn >= 1)) {
    stop("Please use an integer value for Burn (1<=Burn<N); multiple of thin.\n")
  }
  if (!is.logical(Regression)) {
    stop("Please use a logical value for the Regression parameter.\n")
  }

  if (!(PriorMu %in% c("default", "EmpiricalBayes"))) {
    stop("Invalid value for PriorMu. Please use 'default' or 'EmpiricalBayes'")
  }
}

.BASiCS_MCMC_GlobalParams <- function(Data) {
  BioCounts <- counts(Data)
  n <- ncol(Data)
  q.bio <- nrow(Data)

  # If data contains spike-ins
  if (length(altExpNames(Data)) > 1) {
    q <- q.bio + nrow(altExp(Data))
  } else {
    q <- q.bio
  }

  # If Data contains batches
  if (!is.null(colData(Data)$BatchInfo)) {
    # Store the correct number of levels in batch vector
    BatchInfo <- as.factor(colData(Data)$BatchInfo)
    BatchInfo <- factor(BatchInfo, levels = unique(BatchInfo))
    nBatch <- length(unique(BatchInfo))
  } else {
    BatchInfo <- rep(1, times = n)
    nBatch <- 1
  }

  # Parameters associated to the presence of batches
  if (nBatch > 1) {
    BatchDesign <- model.matrix(~ BatchInfo - 1)
  } else {
    # If there are no batches or the user asked to ignore them
    BatchDesign <- matrix(1, nrow = n, ncol = 1)
  }

  list(
    BioCounts = as.matrix(BioCounts), q = q, q.bio = q.bio, n = n,
    nBatch = nBatch, BatchInfo = BatchInfo, BatchDesign = BatchDesign
  )
}

.EmpiricalBayesMu <- function(Data, PriorParam) {

  # Apply scran normalisation
  s <- calculateSumFactors(Data)
  NormCounts <- t(t(counts(Data)) / s)
  aux <- rowMeans(NormCounts)

  # Correction for those genes with total count = 0
  aux <- ifelse(aux == 0, min(aux[aux > 0]), aux)

  PriorParam$mu.mu <- log(aux) - PriorParam$s2.mu / 2

  return(PriorParam)
}

.BASiCS_MCMC_ExtraArgs <- function(Data,
                                   Burn,
                                   GPar,
                                   Regression,
                                   WithSpikes,
                                   PriorMu = c("default", "EmpiricalBayes"),
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
                                     mu.mu = rep(0, times = nrow(Data)),
                                     s2.mu = 0.5,
                                     s2.delta = 0.5,
                                     a.delta = 1,
                                     b.delta = 1,
                                     p.phi = rep(1, times = GPar$n),
                                     a.s = 1,
                                     b.s = 1,
                                     a.theta = 1,
                                     b.theta = 1,
                                     GeneExponent = 1,
                                     CellExponent = 1
                                   ),
                                   eta = 5,
                                   Start = HiddenBASiCS_MCMC_Start(
                                     Data,
                                     PriorParam,
                                     WithSpikes,
                                   ),
                                   mintol_mu = 1e-3,
                                   mintol_delta = 1e-3,
                                   mintol_nu = 1e-5,
                                   mintol_theta = 1e-4) {
  PriorMu <- match.arg(PriorMu)
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

  if (Regression) {
    PriorParam$m <- rep(0, k)
    PriorParam$V <- diag(k)
    PriorParam$a.sigma2 <- 2
    PriorParam$b.sigma2 <- 2
    PriorParam$eta <- eta
  }

  # Redefine prior under the empirical Bayes approach
  if (PriorMu == "EmpiricalBayes") {
    PriorParam <- .EmpiricalBayesMu(Data, PriorParam)
  }

  # Validity checks
  if (
    !(
      # PriorParam$mu.mu > 0 &
      length(PriorParam$mu.mu) == nrow(counts(Data)) &
        PriorParam$s2.mu > 0 & length(PriorParam$s2.mu) == 1 &
        PriorParam$s2.delta > 0 & length(PriorParam$s2.delta) == 1 &
        PriorParam$a.delta > 0 & length(PriorParam$a.delta) == 1 &
        PriorParam$b.delta > 0 & length(PriorParam$b.delta) == 1 &
        all(PriorParam$p.phi > 0) & length(PriorParam$p.phi) == GPar$n &
        PriorParam$a.s > 0 & length(PriorParam$a.s) == 1 &
        PriorParam$b.s > 0 & length(PriorParam$b.s) == 1 &
        PriorParam$a.theta > 0 & length(PriorParam$a.theta) == 1 &
        PriorParam$b.theta > 0 & length(PriorParam$b.theta) == 1 &
        PriorParam$GeneExponent > 0 & length(PriorParam$GeneExponent) == 1 &
        PriorParam$CellExponent > 0 & length(PriorParam$CellExponent) == 1
    )
  ) {
    stop("Invalid prior hyper-parameter values.")
  }
  if (Regression) {
    if (!(length(PriorParam$m) == k &
      ncol(PriorParam$V) == k & nrow(PriorParam$V) == k &
      PriorParam$a.sigma2 > 0 & length(PriorParam$a.sigma2) == 1 &
      PriorParam$b.sigma2 > 0 & length(PriorParam$b.sigma2) == 1)) {
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

  if (!(mintol_mu > 0) | !(mintol_delta > 0) |
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
    ConstrainProp
  )
  ConstrainGene <- NoSpikesParam$ConstrainGene
  NotConstrainGene <- NoSpikesParam$NotConstrainGene
  Constrain <- NoSpikesParam$Constrain
  RefGenes <- NoSpikesParam$RefGenes
  RefGene <- NoSpikesParam$RefGene
  Index <- seq_len(GPar$q.bio) - 1

  list(
    AR = AR, StopAdapt = StopAdapt, StoreChains = StoreChains,
    StoreAdapt = StoreAdapt, StoreDir = StoreDir,
    RunName = RunName, PrintProgress = PrintProgress,
    PriorParam = PriorParam, PriorDeltaNum = PriorDeltaNum,
    PriorDelta = PriorDelta,
    StochasticRef = StochasticRef,
    ConstrainType = ConstrainType, ConstrainProp = ConstrainProp,
    k = k, variance = variance, Start = Start,
    mintol_mu = mintol_mu, mintol_delta = mintol_delta,
    mintol_nu = mintol_nu, mintol_theta = mintol_theta,
    ConstrainGene = ConstrainGene, NotConstrainGene = NotConstrainGene,
    Constrain = Constrain, RefGenes = RefGenes, RefGene = RefGene, Index = Index
  )
}
