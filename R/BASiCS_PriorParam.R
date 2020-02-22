#' Prior parameters for BASiCS_MCMC
#'
#' This is a convenience function to allow partial specification of prior 
#' parameters, and to ensure default parameters are consistent across usage
#' within the package.
#' 
#' \describe{
#'   \item{
#'     \code{mu.mu}
#'   }{
#'     Mean hyper-parameter for the
#'     log-Normal(\code{mu.mu},\code{s2.mu}) prior that is shared by all
#'     gene-specific expression rate parameters \eqn{\mu_i}.
#'     Default: \code{s2.mu = 0}.
#'   }
#'   \item{
#'     \code{s2.mu}
#'   }{
#'     Scale hyper-parameter for the
#'     log-Normal(\code{mu.mu},\code{s2.mu}) prior that is shared by all
#'     gene-specific expression rate parameters \eqn{\mu_i}.
#'     Default: \code{s2.mu = 0.5}.
#'   }
#'   \item{
#'     \code{s2.delta}
#'   }{
#'     Only used when `PriorDelta == 'log-normal'`.
#'     Scale hyper-parameter for the log-Normal(\code{0},\code{s2.delta})
#'     prior that is shared by all gene-specific over-dispersion parameters
#'     \eqn{\delta_i}. Default: \code{s2.delta = 0.5}.
#'   }
#'   \item{
#'     \code{a.delta}
#'   }{
#'     Only used when `PriorDelta == 'gamma'`.
#'     Shape hyper-parameter for the Gamma(\code{a.delta},\code{b.delta})
#'     prior that is shared by all gene-specific biological over-dispersion
#'     parameters \eqn{\delta_i}. Default: \code{a.delta = 1}.
#'   }
#'   \item{
#'     \code{b.delta}
#'   }{
#'     Only used when `PriorDelta == 'gamma'`.
#'     Rate hyper-parameter for the Gamma(\code{a.delta},\code{b.delta})
#'     prior that is shared by all gene-specific biological over-dispersion
#'     hyper-parameters \eqn{\delta_i}. Default: \code{b.delta = 1}.
#'   }
#'   \item{
#'     \code{p.phi}
#'   }{
#'     Dirichlet hyper-parameter for the joint of all
#'     (scaled by \code{n}) cell-specific mRNA content normalising
#'     constants \eqn{\phi_j / n}.
#'     Default: \code{p.phi} \code{= rep(1, n)}.
#'   }
#'   \item{
#'     \code{a.s}
#'   }{
#'     Shape hyper-parameter for the
#'     Gamma(\code{a.s},\code{b.s}) prior that is shared by all
#'     cell-specific capture efficiency normalising constants \eqn{s_j}.
#'     Default: \code{a.s = 1}.
#'   }
#'   \item{
#'     \code{b.s}
#'   }{
#'     Rate hyper-parameter for the Gamma(\code{a.s},
#'     \code{b.s}) prior that is shared by all cell-specific capture
#'     efficiency normalising constants \eqn{s_j}.
#'     Default: \code{b.s = 1}.
#'   }
#'   \item{
#'     \code{a.theta}
#'   }{
#'     Shape hyper-parameter for the
#'     Gamma(\code{a.theta},\code{b.theta}) prior for technical noise
#'     parameter \eqn{\theta}. Default: \code{a.theta = 1}.
#'   }
#'   \item{
#'     \code{b.theta}
#'   }{
#'     Rate hyper-parameter for the
#'     Gamma(\code{a.theta},\code{b.theta}) prior for technical noise
#'     parameter \eqn{\theta}. Default: \code{b.theta = 1}.
#'   }
#'   \item{
#'     \code{eta}
#'   }{
#'     Only used when \code{Regression = TRUE}. \code{eta}
#'     specifies the degress of freedom for the residual term.
#'     Default: \code{eta = 5}.
#'   }
#'   \item{
#'     \code{k}
#'   }{
#'     Only used when \code{Regression = TRUE}. \code{k} specifies
#'     the number of regression Gaussian Radial Basis Functions (GRBF) used
#'     within the correlated prior adopted for gene-specific over-dispersion
#'     and mean expression parameters. Default: \code{k = 12}.
#'   }
#' }
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
#' during burn-in?
#' @param RBFLocations RBFLocations of GRBFs in units of \code{log(mu)}
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
    FixLocations = !is.null(RBFLocations),
    RBFLocations = NULL,
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
    RBFLocations = RBFLocations,
    m = m,
    V = V,
    k = k,
    a.sigma2 = a.sigma2,
    b.sigma2 = b.sigma2,
    eta = eta
  )
}
