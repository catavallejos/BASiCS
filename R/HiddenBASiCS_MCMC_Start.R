# Used in BASiCS_MCMC
HiddenBASiCS_MCMC_Start <- function(Data,
                                    PriorParam,
                                    WithSpikes,
                                    ls.mu0 = -4,
                                    ls.delta0 = -2,
                                    ls.phi0 = 11,
                                    ls.nu0 = -10,
                                    ls.theta0 = -4)
{
  if (!is(Data, "SingleCellExperiment")) {
    stop("'Data' is not a SingleCellExperiment class object.")
  }

  CountsBio <- counts(Data)
  n <- ncol(Data)
  q.bio <- nrow(CountsBio)

  # Extract spike-in genes
  if (WithSpikes) { CountsTech <- assay(altExp(Data)) }

  # Initialize normalization as the 'scran' estimates
  suppressWarnings(size_scran <- scran::computeSumFactors(CountsBio))
  # Fix for cases in which 'scran' normalisation has invalid output
  if ((min(size_scran) <= 0) | (sum(is.na(size_scran)) > 0)) {
    message("-------------------------------------------------------------\n",
            "There was an issue when applying `scran` normalization  \n",
            "`positive = TRUE` has been added to `computeSumFactors` call \n",
            "Please consider a more stringent quality control criteria. \n",
            "-------------------------------------------------------------\n")
    suppressWarnings(size_scran <- scran::computeSumFactors(CountsBio,
                                                            positive = TRUE))
  }

  if (WithSpikes) {
    # Initialize s as the empirical capture efficiency rates
    s0 <- matrixStats::colSums2(CountsTech) /
      sum(metadata(Data)$SpikeInput[,2])
    nu0 <- s0
    phi0 <- size_scran / s0
    phi0 <- n * phi0 / sum(phi0)

    # Initialize mu using average 'normalised counts' across cells
    # and true input values for spike-in genes
    nCountsBio <- t( t(CountsBio) / (phi0 * s0) )
    meansBio <- rowMeans(nCountsBio)
    # +1 to avoid zeros as starting values
    #mu0 <- c(meansBio + 1, metadata(Data)$SpikeInput)
    mu0 <- meansBio + 1
  }
  else {
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

    # Starting value for delta
    # Defined by the CV for high- and mid-expressed genes
    # This is motivated by equation (2) in Vallejos et al (2016)
    varsBio <- matrixStats::rowVars(nCountsBio)
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
    out <- list(mu0 = mu0, delta0 = delta0,
                phi0 = phi0, s0 = s0,
                nu0 = nu0, theta0 = theta0,
                ls.mu0 = ls.mu0, ls.delta0 = ls.delta0,
                ls.phi0 = ls.phi0, ls.nu0 = ls.nu0, ls.theta0 = ls.theta0)

    # Starting values for regression-related parameters
    if (!is.null(PriorParam$eta)) {
      out$beta0 <- mvrnorm(1, PriorParam$m, PriorParam$V)
      out$sigma20 <- rgamma(1, PriorParam$a.sigma2, PriorParam$b.sigma2)
      out$lambda0 <- rgamma(q.bio, shape = PriorParam$eta/2,
                        rate = PriorParam$eta/2)
    }

    return(out)
}

