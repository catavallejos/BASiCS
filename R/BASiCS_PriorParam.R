#' Prior parameters for BASiCS_MCMC
#'
#' This is a convenience function to allow partial specification of prior 
#' parameters, and to ensure default parameters are consistent across usage
#' within the package.
#' 
#' @param Data \linkS4class{SingleCellExperiment} object (required).
#' @param k Number of regression terms, including k - 2 Gaussian 
#' radial basis functions (GRBFs).
#' @param mu.mu,s2.mu Mean and variance parameters for lognormal prior on mu.
#' @param s2.delta Variance parameter for lognormal prior on delta when 
#' \code{PriorDelta="lognormal"}.
#' @param a.delta,b.delta Parameters for gamma prior on delta when 
#' \code{PriorDelta="gamma"}.
#' @param p.phi Parameter for dirichlet prior on phi.
#' @param a.s,b.s Parameters for gamma prior on s.
#' @param a.theta,b.theta Parameters for gamma prior on theta.
#' @param RBFMinMax Should GRBFs be placed at the minimum and maximum of 
#' \code{log(mu)}?
#' @param FixLocations Should RBFLocations be fixed throughout MCMC, or adaptive
#' during burn-in? By default this is \code{FALSE}, but it is set to \code{TRUE}
#' if \code{RBFLocations} or \code{MinGenesPerRBF} are specified.
#' @param RBFLocations Numeric vector specifying locations of GRBFs in units 
#' of \code{log(mu)}.
#' @param MinGenesPerRBF Numeric scalar specifying the minimum number of genes
#' for GRBFs to be retained. If fewer than \code{MinGenesPerRBF} genes have 
#' values of \code{log(mu)} within the range of an RBF, it is removed. 
#' The range covered by each RBF is defined as centre of the RBF plus or minus
#' half the distance between RBFs.
#' @param variance Variance of the GRBFs.
#' @param m,V Mean and (co)variance priors for regression coefficients.
#' @param a.sigma2,b.sigma2 Priors for inverse gamma prior on regression scale.
#' @param eta Degrees of freedom for t distribution of regresion errors.
#' @param PriorMu Indicates if the original prior (\code{PriorMu = 'default'})
#' or an empirical Bayes approach (\code{PriorMu = 'EmpiricalBayes'}) will be 
#' assigned to gene-specific mean expression parameters.
#' @param PriorDelta Scalar character specifying the prior type to use for 
#' delta overdispersion parameter. Options are "log-normal" (recommended)
#' and "gamma" (used in Vallejos et al. (2015)).
#' @param StochasticRef Logical scalar specifying whether the reference gene
#' for the no-spikes version should be chosen randomly at MCMC iterations.
#' @param ConstrainProp Proportion of genes to be considered as reference genes
#' if \code{StochasticRef=TRUE}.
#' @param GeneExponent,CellExponent Exponents for gene and cell-specific 
#' parameters. These should not be outside of divide and conquer MCMC 
#' applications.
#' 
#' @return A list containing the prior hyper-parameters that are required to
#' run the algoritm implemented in \code{\link[BASiCS]{BASiCS_MCMC}}.
#'
#' @examples
#' 
#' BASiCS_PriorParam(makeExampleBASiCS_Data(), k = 12)
#' 
#' 
#' @export
BASiCS_PriorParam <- function(
    Data,
    k = 12,
    mu.mu = if (PriorMu == "default") rep(0, nrow(Data))
      else .EmpiricalBayesMu(Data, s2.mu, with_spikes = WithSpikes),
    s2.mu = 0.5,
    s2.delta = 0.5,
    a.delta = 1, 
    b.delta = 1,
    p.phi = rep(1, times = ncol(Data)),
    a.s = 1,
    b.s = 1,
    a.theta = 1,
    b.theta = 1,
    RBFMinMax = TRUE,
    FixLocations = !is.null(RBFLocations) | !is.na(MinGenesPerRBF),
    RBFLocations = NULL,
    MinGenesPerRBF = NA,
    variance = 1.2,
    m = numeric(k),
    V = diag(k),
    a.sigma2 = 2,
    b.sigma2 = 2,
    eta = 5,
    PriorMu = c('default', 'EmpiricalBayes'),
    PriorDelta = c("log-normal", "gamma"),
    StochasticRef = TRUE,
    ConstrainProp = 0.2,
    WithSpikes = FALSE,
    GeneExponent = 1,
    CellExponent = 1
  ) {
  
 # Removed as I do not expect any users to select this option  
 # if (!missing(ConstrainType)) {
 #   warning("Use of ConstrainType is deprecated.")
 # }
  PriorMu <- match.arg(PriorMu)
  PriorDelta <- match.arg(PriorDelta)
  n <- ncol(Data)
  q <- nrow(Data)
  list(
    mu.mu = mu.mu,
    s2.mu = s2.mu,
    s2.delta = s2.delta,
    a.delta = a.delta,
    b.delta = b.delta,
    p.phi = p.phi,
    GeneExponent = GeneExponent,
    CellExponent = CellExponent,
    PriorMu = PriorMu,
    PriorDelta = PriorDelta,
    StochasticRef = StochasticRef,
    ConstrainProp = ConstrainProp,
    MinGenesPerRBF = MinGenesPerRBF,
    a.s = a.s,
    b.s = b.s,
    a.theta = a.theta,
    b.theta = b.theta,
    RBFMinMax = RBFMinMax,
    FixLocations = FixLocations,
    RBFLocations = RBFLocations,
    variance = variance,
    m = m,
    V = V,
    k = k,
    a.sigma2 = a.sigma2,
    b.sigma2 = b.sigma2,
    eta = eta
  )
}
