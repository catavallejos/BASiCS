# Used in BASiCS_MCMC
HiddenBASiCS_MCMC_Start <- function(Data, 
                                    k = NULL, variance = NULL, eta = NULL,
                                    ...) 
{
  if (!is(Data, "SingleCellExperiment")) 
    stop("'Data' is not a SingleCellExperiment class object.")
  
  # Number of instrinsic genes
  q <- length(isSpike(Data))
  q.bio <- sum(!isSpike(Data))
  # Number of cells
  n <- dim(assay(Data))[2]
  
  # Separating spike-ins from the rest of genes
  CountsBio <- as.matrix(assay(Data)[!isSpike(Data), , drop = FALSE]) 
  CountsTech <- as.matrix(assay(Data)[isSpike(Data), , drop = FALSE])
  
  # Initialize normalization as the 'scran' estimates
  # Change default pool sizes for small sample size
  sizes.aux <- c(20, 40, 60, 80, 100)
  if (n < 200) { sizes.aux <- c(20, 40, 60, 80) }
  if (n < 160) { sizes.aux <- c(20, 40, 60) }
  if (n < 120) { sizes.aux <- c(20, 40) }
  if (n < 80) { sizes.aux <- c(20) }
  if (n < 40) { sizes.aux <- c(10) }
  
  size_scran <- scran::computeSumFactors(CountsBio, sizes = sizes.aux)
  
  if (length(metadata(Data)$SpikeInput) > 1) 
  {
    # Initialize s as the empirical capture efficiency rates
    s0 <- matrixStats::colSums2(CountsTech) / 
      sum(metadata(Data)$SpikeInput)
    nu0 <- s0
    phi0 <- size_scran / s0
    phi0 <- n * phi0 / sum(phi0)
    
    # Initialize mu using average 'normalised counts' across cells 
    # and true input values for spike-in genes
    nCountsBio <- t( t(CountsBio) / (phi0 * s0) )
    meansBio <- rowMeans(nCountsBio)
    # +1 to avoid zeros as starting values
    mu0 <- c(meansBio + 1, metadata(Data)$SpikeInput)  
  } 
  else 
  {
    s0 <- size_scran
    s0 <- n * s0 / sum(s0)
    for (B in unique(metadata(Data)$BatchInfo)) 
    {
      aux <- metadata(Data)$BatchInfo == B
      s0[aux] <- sum(aux) * s0[aux]/sum(s0[aux])
    }
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
    cv2Bio <- varsBio/(meansBio)^2
    delta0 <- rgamma(q.bio, 1, 1) + 1
    Aux <- which(meansBio > stats::quantile(meansBio, 0.1))
    delta0[Aux] <- cv2Bio[Aux]
    # 1e-3 added to be coherent with tolerance used within MCMC sampler
    delta0 <- delta0 + 0.001
    
    # Random stating value for theta (within typically observed range)
    theta0 <- runif(1, min = 0.2, max = 1)
    
    # If given, load default values for adaptive proposal variances
    args <- list(...)
    ls.mu0 <- ifelse("ls.mu0" %in% names(args), args$ls.mu0, -4)
    ls.delta0 <- ifelse("ls.delta0" %in% names(args), args$ls.delta0, -2)
    ls.phi0 <- ifelse("ls.phi0" %in% names(args), args$ls.phi0, 11)
    ls.nu0 <- ifelse("ls.nu0" %in% names(args), args$ls.nu0, -10)
    ls.theta0 <- ifelse("ls.theta0" %in% names(args), args$ls.theta0, -4)
    
    # Starting values for regression approach
    if(!is.null(k)){
      m0 = rep(0, k); V0 = diag(k); sigma2.a0 = 2; sigma2.b0 = 2
      beta0 = mvrnorm(1,m0,V0); sigma20 = rgamma(1,sigma2.a0,sigma2.b0)
      reg.nu0 = eta
      lambda0 = rgamma(q.bio,shape=reg.nu0/2,rate=reg.nu0/2)
    }
    
    # Starting values for the proposal variances 
    ls.mu0 <- rep(ls.mu0, q.bio)
    ls.delta0 <- rep(ls.delta0, q.bio)
    ls.phi0 <- ifelse(n < 200, pmax(2 * log(n), ls.phi0), 11) 
    ls.nu0 <- pmax(2 * log(0.02 * abs(log(nu0))), ls.nu0)
    ls.theta0 <- pmax(2 * log(0.02 * abs(log(theta0))), ls.theta0)
    
    if(!is.null(k)){
      list(mu0 = mu0, delta0 = delta0, 
           phi0 = phi0, s0 = s0, 
           nu0 = nu0, theta0 = theta0, 
           ls.mu0 = ls.mu0, ls.delta0 = ls.delta0, 
           ls.phi0 = ls.phi0, ls.nu0 = ls.nu0, ls.theta0 = ls.theta0,
           m0 = m0, V0 = V0, sigma2.a0 = sigma2.a0, sigma2.b0 = sigma2.b0,
           beta0 = beta0, sigma20 = sigma20,
           lambda0 = lambda0)
    }
    else{
      list(mu0 = mu0, delta0 = delta0, 
           phi0 = phi0, s0 = s0, 
           nu0 = nu0, theta0 = theta0, 
           ls.mu0 = ls.mu0, ls.delta0 = ls.delta0, 
           ls.phi0 = ls.phi0, ls.nu0 = ls.nu0, ls.theta0 = ls.theta0)
    }
}

