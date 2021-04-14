.BASiCS_MCMC_InputCheck <- function(Data,
                                    N,
                                    Thin,
                                    Burn,
                                    Regression,
                                    WithSpikes,
                                    Threads,
                                    NSubsets) {

  if (!is(Data, "SingleCellExperiment")) {
    stop("'Data' is not a SingleCellExperiment class object.\n")
  }
  if (!is.numeric(Threads) &&
      length(Threads) == 1 &&
      round(Threads) == Threads &&
      Threads > 0) {
    stop("Threads must be a positive integer.")
  }

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
      "altExp '", altExpNames(Data),"' is assumed to contain spike-in genes.\n",
      "see help(altExp) for details. \n"
    )
    # If SpikeInput slot is missing and WithSpikes == TRUE
    if (is.null(rowData(altExp(Data)))) {
      stop(
        "'altExp(Data)' does not contain 'rowData' \n"
      )  
    }
  }
  
  # If SpikeInput slot is missing and WithSpikes == TRUE
  if (length(altExpNames(Data)) > 0 &&
     WithSpikes &&
     is.null(rowData(altExp(Data)))) {
    stop(
      "'altExp(Data)' does not contain 'rowData' \n"
    )
  }
  
  # If isSpike slot is missing and WithSpikes == TRUE
  if ((length(altExpNames(Data)) == 0)  & WithSpikes) {
    stop(
      "'Data' does not contain information about spike-in genes \n", 
      "Please include this information using 'altExp' \n",
      "or set 'WithSpikes = FALSE' \n."
    )
  }
    
  # If BatchInfo slot is missing and WithSpikes == FALSE
  if (!WithSpikes & is.null(colData(Data)$BatchInfo)) {
    warning(
      "'Data' should contain a BatchInfo vector when 'WithSpikes = FALSE'. \n", 
      "Please assign the batch information to: \n
      'colData(Data)$BatchInfo = BatchInfo'. \n"
    )
  }

  # Checking how counts are stored
  if (!("counts" %in% assayNames(Data)))
    stop(
      "'Data' does not contain a 'counts' slot. \n",
      "Please make sure to include the raw data in the \n", 
      "SingleCellExperiment object under the name 'counts'."
    )
  
  # Further checks on input
  errors <- .ChecksBASiCS_Data(Data, WithSpikes)
  
  if (length(errors) > 0) {
    stop(errors) 
  }
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
    stop("Burn should be an integer; 1<=Burn<N; multiple of thin.\n")
  }
  if (!is.logical(Regression)) {
    stop("Please use a logical value for the Regression parameter.\n")  
  }
  if (!is.numeric(NSubsets) ||
      length(NSubsets) > 1 ||
      NSubsets < 1 ||
      round(NSubsets) != NSubsets) {
    stop("Invalid value for NSubsets; should be a length 1 positive integer")
  }
}

.BASiCS_MCMC_Start <- function(Data,
                               PriorParam,
                               WithSpikes,
                               Regression,
                               ls.mu0 = -4,
                               ls.delta0 = -2,
                               ls.phi0 = 11,
                               ls.nu0 = -10,
                               ls.theta0 = -4) {
  
  if (!is(Data, "SingleCellExperiment")) {
    stop("'Data' is not a SingleCellExperiment class object.")
  }

  CountsBio <- counts(Data)
  n <- ncol(Data)
  q.bio <- nrow(CountsBio)

  # Extract spike-in genes
  if (WithSpikes) {
    CountsTech <- assay(altExp(Data))
  }

  # Initialize normalization as the 'scran' estimates
  suppressWarnings(size_scran <- scran::calculateSumFactors(CountsBio))
  # Fix for cases in which 'scran' normalisation has invalid output
  if ((min(size_scran) <= 0) | (sum(is.na(size_scran)) > 0)) {
    message(
      "-------------------------------------------------------------\n",
      "There was an issue when applying `scran` normalization  \n",
      "Running scran::cleanSizeFactors on size factors to ensure positivity.\n",
      "Please consider a more stringent quality control criteria. \n",
      "-------------------------------------------------------------\n"
    )
    size_scran <- scran::cleanSizeFactors(
      size_scran,
      num.detected = colSums(CountsBio != 0)
    )
  }

  if (WithSpikes) {
    # Initialize s as the empirical capture efficiency rates
    s0 <- Matrix::colSums(CountsTech) /
      sum(rowData(altExp(Data))[, 2])
    nu0 <- s0
    phi0 <- size_scran / s0
    phi0 <- n * phi0 / sum(phi0)

    # Initialize mu using average 'normalised counts' across cells
    # and true input values for spike-in genes
    nCountsBio <- t(t(CountsBio) / (phi0 * s0))
    meansBio <- rowMeans(nCountsBio)
    # +1 to avoid zeros as starting values
    mu0 <- meansBio + 1
  } else {
    s0 <- size_scran
    nu0 <- s0
    phi0 <- NULL

    # Initialize mu using average 'normalised counts' across cells
    nCountsBio <- t( t(CountsBio) / s0 )
    meansBio <- rowMeans(nCountsBio)
    # +1 to avoid zeros as starting values
    meansBio <- ifelse(meansBio == 0, meansBio + 1, meansBio)
    mu0 <- meansBio
  }
  
  # If EB prior, replace mu0 by EB estimate
  if (length(unique(PriorParam$mu.mu)) > 1) {
    mu0 <- .EmpiricalBayesMu(
      Data,
      s2_mu = PriorParam$s2.mu,
      with_spikes = WithSpikes,
      log_scale = FALSE
    )
  }

  # Starting value for delta
  # Defined by the CV for high- and mid-expressed genes
  # This is motivated by equation (2) in Vallejos et al (2016)
  varsBio <- apply(nCountsBio, 1, var)
  cv2Bio <- varsBio / (meansBio)^2
  delta0 <- rgamma(q.bio, 1, 1) + 1
  Aux <- which(meansBio > stats::quantile(meansBio, 0.1))
  delta0[Aux] <- cv2Bio[Aux]
  # 1e-3 added to be coherent with tolerance used within MCMC sampler
  delta0 <- delta0 + 0.001

  # Random stating value for theta (within typically observed range)
  theta0 <- runif(1, min = 0.2, max = 1)

  # Starting values for the proposal variances
  ls.mu0 <- rep(ls.mu0, q.bio)
  ls.delta0 <- rep(ls.delta0, q.bio)
  ls.phi0 <- ifelse(n < 200, pmax(2 * log(n), ls.phi0), 11)
  ls.nu0 <- pmax(2 * log(0.02 * abs(log(nu0))), ls.nu0)
  ls.theta0 <- pmax(2 * log(0.02 * abs(log(theta0))), ls.theta0)
  # Convert to numeric values
  ls.phi0 <- as.numeric(ls.phi0)
  ls.theta0 <- as.numeric(ls.theta0)

  # Output list
  out <- list(
    mu0 = mu0,
    delta0 = delta0,
    phi0 = phi0,
    s0 = s0,
    nu0 = nu0,
    theta0 = theta0,
    ls.mu0 = ls.mu0,
    ls.delta0 = ls.delta0,
    ls.phi0 = ls.phi0,
    ls.nu0 = ls.nu0,
    ls.theta0 = ls.theta0
  )

  # Starting values for regression-related parameters
  if (Regression) {
    out$beta0 <- mvrnorm(1, PriorParam$m, PriorParam$V)
    out$sigma20 <- rgamma(1, PriorParam$a.sigma2, PriorParam$b.sigma2)
    out$lambda0 <- rgamma(
      q.bio,
      shape = PriorParam$eta / 2,
      rate = PriorParam$eta / 2
    )
  }

  return(out)
}

.BASiCS_MCMC_GlobalParams <- function(Data) {
  
  BioCounts <- counts(Data)
  n <- ncol(Data)
  q.bio <- nrow(Data)
  
  # If data contains spike-ins
  if (length(altExpNames(Data)) > 1){
    q <- q.bio + nrow(altExp(Data))
  } else { q <- q.bio }
  
  # If Data contains batches
  if(!is.null(colData(Data)$BatchInfo)){
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
    BatchDesign <- model.matrix(~BatchInfo - 1)  
  } else { 
    # If there are no batches or the user asked to ignore them
    BatchDesign <- matrix(1, nrow = n, ncol = 1) 
  }
  
  list(
    BioCounts = as.matrix(BioCounts),
    q = q,
    q.bio = q.bio,
    n = n,
    nBatch = nBatch,
    BatchInfo = BatchInfo,
    BatchDesign = BatchDesign
  )
}

.EmpiricalBayesMu <- function(Data, s2_mu, with_spikes, log_scale = TRUE) {

  # Apply scran normalisation
  s <- calculateSumFactors(Data)
  NormCounts <- t(t(counts(Data)) / s)
  aux <- rowMeans(NormCounts)

  # Correction for those genes with total count = 0
  aux <- ifelse(aux == 0, min(aux[aux > 0]), aux)
  
  # Extra scaling if the data has spike-ins
  if (with_spikes) {
    # overall capture for spike-ins across a cell
    s0 <- Matrix::colSums(assay(altExp(Data))) /
      sum(rowData(altExp(Data))[, 2])  
    # arbitrary scaling set to match the mode of the distribution of s0
    s0_dens <- density(s0)
    myscale <- s0_dens$x[s0_dens$y == max(s0_dens$y)]
    aux <- aux / myscale
  }

  if (log_scale) {
    log(aux) - s2_mu / 2
  } else {
    aux * exp(-s2_mu / 2)
  }
}

.BASiCS_MCMC_ExtraArgs <- function(Data,
                                   Burn,
                                   GPar,
                                   Regression,
                                   WithSpikes,
                                   AR = 0.44,
                                   StopAdapt = Burn,
                                   StoreChains = FALSE,
                                   StoreAdapt = FALSE,
                                   StoreDir = getwd(),
                                   RunName = "",
                                   PrintProgress = TRUE,
                                   PriorParam = BASiCS_PriorParam(Data, k = 12),
                                   Start = .BASiCS_MCMC_Start(
                                     Data = Data,
                                     PriorParam = PriorParam,
                                     WithSpikes = WithSpikes,
                                     Regression = Regression
                                   ),
                                   mintol_mu = 1e-5,
                                   mintol_delta = 1e-3,
                                   mintol_nu = 1e-5,
                                   mintol_theta = 1e-4) {

  .stop_k(PriorParam$k)
  lm <- log(Start$mu0)
  if (is.null(PriorParam$RBFLocations)) {
    RBFLocations <- .estimateRBFLocations(
      lm,
      PriorParam$k,
      PriorParam$RBFMinMax
    )
  } else {
    RBFLocations <- PriorParam$RBFLocations
  }
  if (!is.na(PriorParam$MinGenesPerRBF)) {
    d <- (RBFLocations[[2]] - RBFLocations[[1]]) / 2
    retain <- vapply(
      RBFLocations,
      function(RBFLocation) {
        ind_within <- lm > RBFLocation - d & lm < RBFLocation + d
        sum(ind_within) > PriorParam$MinGenesPerRBF
      },
      logical(1)
    )
    RBFLocations <- RBFLocations[retain]
    ## Intercept & linear term
    retain_prior <- c(TRUE, TRUE, retain)
    PriorParam$m <- PriorParam$m[retain_prior]
    PriorParam$V <- PriorParam$V[retain_prior, retain_prior]
    Start$beta0 <- Start$beta0[retain_prior]
  }
  PriorParam$RBFLocations <- RBFLocations
  PriorParam$k <- length(RBFLocations) + 2

  if (!Regression) {
    message(
      "-------------------------------------------------------------\n",
      "NOTE: default choice PriorDelta = 'log-normal'  (recommended value). \n",
      "Vallejos et al (2015) used a 'gamma' prior instead.\n",
      "-------------------------------------------------------------\n"
    )
  }

  if (Regression) {
    .stop_k(PriorParam$k)
  }

  if (is.null(PriorParam$mu.mu)) {
    PriorParam$mu.mu <- if (PriorParam$PriorMu == "default") rep(0, nrow(Data))
      else .EmpiricalBayesMu(Data, PriorParam$s2.mu, with_spikes = WithSpikes)
  }

  # Validity checks
  assertthat::assert_that(
    length(PriorParam$mu.mu) == nrow(Data),
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
    PriorParam$GeneExponent > 0,
    length(PriorParam$GeneExponent) == 1,
    PriorParam$CellExponent > 0,
    length(PriorParam$CellExponent) == 1,
    StopAdapt > 0,
    length(StoreChains) == 1,
    length(StoreAdapt) == 1,
    dir.exists(StoreDir),
    mintol_mu > 0,
    mintol_delta > 0,
    mintol_nu > 0,
    mintol_theta > 0,
    is.logical(PriorParam$RBFMinMax),
    length(PriorParam$RBFMinMax) == 1,
    is.logical(PriorParam$FixLocations),
    is.na(PriorParam$MinGenesPerRBF) || (
      length(PriorParam$MinGenesPerRBF) &
      is.numeric(PriorParam$MinGenesPerRBF)
    ),
    length(PriorParam$FixLocations) == 1
  )
  assertthat::assert_that(
    length(Start$mu0) == GPar$q.bio,
    length(Start$delta0) == GPar$q.bio,
    is.null(Start$phi0) || length(Start$phi) == GPar$n,
    length(Start$nu0) == GPar$n,
    length(Start$s0) == GPar$n
    # , length(Start$theta0) == GPar$nBatch
  )

  if (Regression) {
    assertthat::assert_that(
      length(PriorParam$m) == PriorParam$k,
      ncol(PriorParam$V) == PriorParam$k,
      nrow(PriorParam$V) == PriorParam$k,
      PriorParam$a.sigma2 > 0,
      length(PriorParam$a.sigma2) == 1,
      PriorParam$b.sigma2 > 0,
      length(PriorParam$b.sigma2) == 1
    )
  }
  if (!(AR > 0 & AR < 1 & length(AR) == 1)) {
    stop("Invalid AR value. Recommended value: AR = 0.44.")
  }

  ## Definition of parameters that are specific to the no-spikes case
  ## if not, empty list because list()[["foo"]] is null
  if (!WithSpikes) {  
    NoSpikesParam <- .BASiCS_MCMC_NoSpikesParam(
      GPar$BioCounts,
      PriorParam$StochasticRef,
      GPar$q.bio,
      Start$mu0,
      PriorParam$PriorDelta,
      PriorParam$ConstrainProp
    )
  } else {
    NoSpikesParam <- list()
  }
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
    Index = Index
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


.BASiCS_MCMC_NoSpikesParam <- function(Counts,
                                       StochasticRef,
                                       q.bio,
                                       mu0,
                                       PriorDelta,
                                       ConstrainProp) {

  if (PriorDelta == "gamma") {
    stop("PriorDelta = 'gamma' is not supported for the no-spikes case")
  }

  ConstrainGene <- which(matrixStats::rowSums2(Counts) >= 1) - 1
  NotConstrainGene <- which(matrixStats::rowSums2(Counts) < 1) - 1
  NonZero <- ConstrainGene + 1
  Constrain <- mean(log(mu0[NonZero]))

  # Whether or not a stochatic reference is used
  # If stochastic, range of possible reference values only includes
  # the nearest 10% genes located around the constrain

  if (StochasticRef) {
    aux.ref <- cbind(ConstrainGene, abs(log(mu0[ConstrainGene+1]) - Constrain))
    aux.ref <- aux.ref[order(aux.ref[, 2]), ]
    # In total ConstrainProp*100% of genes to be used as reference candidates
    CandidateRef <- round(ConstrainProp * q.bio)
    # Fix for the code to run on the synthetic small dataset
    # generaed by makeExample_BASiCS function (less than 200 genes)
    if (length(ConstrainGene) > CandidateRef) {
      RefGenes <- aux.ref[seq_len(CandidateRef), 1]
    } else {
      RefGenes <- aux.ref[, 1]
    }
    RefGene <- RefGenes[1]
  } else {
    aux.ref <- which(abs(log(mu0[ConstrainGene+1]) - Constrain) ==
                       min(abs(log(mu0[ConstrainGene+1]) - Constrain)))[1]
    RefGene <- ConstrainGene[aux.ref]
    RefGenes <- RefGene
  }

  list(
    ConstrainGene = ConstrainGene,
    NotConstrainGene = NotConstrainGene,
    Constrain = Constrain,
    RefGenes = RefGenes,
    RefGene = RefGene
  )
}
