#' Prior parameters for BASiCS_MCMC
#'
#' This is a convenience function to allow partial specification of prior 
#' parameters, and to ensure default parameters are consistent across usage
#' within the package.
#' 
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
#' @param FixLocations Should locations be fixed throughout MCMC, or adaptive
#' during burn-in?
#' @param locations Locations of GRBFs in units of \code{log(mu)}
#' @param m,V Mean and (co)variance priors for regression coefficients.
#' @param a.sigma2,b.sigma2 Priors for inverse gamma prior on regression scale.
#' @param eta Degrees of freedom for t distribution of regresion errors.
#' @param GeneExponent,CellExponent Exponents for gene and cell-specific 
#' parameters. These should not be outside of divide and conquer MCMC 
#' applications.
#'
#' @examples
#' 
#' BASiCS_PriorParam(makeExampleBASiCS_Data(), 12)
#' 
#' @export
BASiCS_PriorParam <- function(
    Data,
    k = 12,
    mu.mu = 0,
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
    FixLocations = FALSE,
    locations = numeric(k - 2),
    m = numeric(k),
    V = diag(k),
    a.sigma2 = 2,
    b.sigma2 = 2,
    eta = 5,
    GeneExponent = 1,
    CellExponent = 1
  ) {

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
    a.s = a.s,
    b.s = b.s,
    a.theta = a.theta,
    b.theta = b.theta,
    RBFMinMax = RBFMinMax,
    FixLocations = FixLocations,
    locations = locations,
    m = m,
    V = V,
    a.sigma2 = a.sigma2,
    b.sigma2 = b.sigma2,
    eta = eta
  )
}
