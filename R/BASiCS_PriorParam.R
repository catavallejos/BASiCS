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
        else .EmpiricalBayesMu(Data, s2.mu),
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


## Old documentation
# \describe{
#   \item{
#     \code{mu.mu}
#   }{
#     Mean hyper-parameter for the
#     log-Normal(\code{mu.mu},\code{s2.mu}) prior that is shared by all
#     gene-specific expression rate parameters \eqn{\mu_i}.
#     Default: \code{s2.mu = 0}.
#   }
#   \item{
#     \code{s2.mu}
#   }{
#     Scale hyper-parameter for the
#     log-Normal(\code{mu.mu},\code{s2.mu}) prior that is shared by all
#     gene-specific expression rate parameters \eqn{\mu_i}.
#     Default: \code{s2.mu = 0.5}.
#   }
#   \item{
#     \code{s2.delta}
#   }{
#     Only used when `PriorDelta == 'log-normal'`.
#     Scale hyper-parameter for the log-Normal(\code{0},\code{s2.delta})
#     prior that is shared by all gene-specific over-dispersion parameters
#     \eqn{\delta_i}. Default: \code{s2.delta = 0.5}.
#   }
#   \item{
#     \code{a.delta}
#   }{
#     Only used when `PriorDelta == 'gamma'`.
#     Shape hyper-parameter for the Gamma(\code{a.delta},\code{b.delta})
#     prior that is shared by all gene-specific biological over-dispersion
#     parameters \eqn{\delta_i}. Default: \code{a.delta = 1}.
#   }
#   \item{
#     \code{b.delta}
#   }{
#     Only used when `PriorDelta == 'gamma'`.
#     Rate hyper-parameter for the Gamma(\code{a.delta},\code{b.delta})
#     prior that is shared by all gene-specific biological over-dispersion
#     hyper-parameters \eqn{\delta_i}. Default: \code{b.delta = 1}.
#   }
#   \item{
#     \code{p.phi}
#   }{
#     Dirichlet hyper-parameter for the joint of all
#     (scaled by \code{n}) cell-specific mRNA content normalising
#     constants \eqn{\phi_j / n}.
#     Default: \code{p.phi} \code{= rep(1, n)}.
#   }
#   \item{
#     \code{a.s}
#   }{
#     Shape hyper-parameter for the
#     Gamma(\code{a.s},\code{b.s}) prior that is shared by all
#     cell-specific capture efficiency normalising constants \eqn{s_j}.
#     Default: \code{a.s = 1}.
#   }
#   \item{
#     \code{b.s}
#   }{
#     Rate hyper-parameter for the Gamma(\code{a.s},
#     \code{b.s}) prior that is shared by all cell-specific capture
#     efficiency normalising constants \eqn{s_j}.
#     Default: \code{b.s = 1}.
#   }
#   \item{
#     \code{a.theta}
#   }{
#     Shape hyper-parameter for the
#     Gamma(\code{a.theta},\code{b.theta}) prior for technical noise
#     parameter \eqn{\theta}. Default: \code{a.theta = 1}.
#   }
#   \item{
#     \code{b.theta}
#   }{
#     Rate hyper-parameter for the
#     Gamma(\code{a.theta},\code{b.theta}) prior for technical noise
#     parameter \eqn{\theta}. Default: \code{b.theta = 1}.
#   }
#   \item{
#     \code{eta}
#   }{
#     Only used when \code{Regression = TRUE}. \code{eta}
#     specifies the degress of freedom for the residual term.
#     Default: \code{eta = 5}.
#   }
#   \item{
#     \code{k}
#   }{
#     Only used when \code{Regression = TRUE}. \code{k} specifies
#     the number of regression Gaussian Radial Basis Functions (GRBF) used
#     within the correlated prior adopted for gene-specific over-dispersion
#     and mean expression parameters. Default: \code{k = 12}.
#   }
# }
