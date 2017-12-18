
#' @title Creates a BASiCS_Chain object from pre-computed MCMC chains
#'
#' @description \code{BASiCS_Chain} creates a \code{\linkS4class{BASiCS_Chain}} 
#' object from pre-computed MCMC chains.
#'
#' @param parameters List of matrices containing MCMC chains for each model 
#' parameter. 
#' \describe{
#' \item{mu}{MCMC chain for gene-specific mean expression parameters 
#' \eqn{\mu_i}, biological genes only 
#' (matrix with \code{q.bio} columns, all elements must be positive numbers)}
#' \item{delta}{MCMC chain for gene-specific biological over-dispersion 
#' parameters \eqn{\delta_i}, biological genes only 
#' (matrix with \code{q.bio} columns, all elements must be positive numbers)}
#' \item{phi}{MCMC chain for cell-specific mRNA content normalisation parameters 
#' \eqn{\phi_j} (matrix with \code{n} columns, all elements must be positive 
#' numbers and the sum of its elements must be equal to \code{n})}.
#' This parameter is only used when spike-in genes are available.
#' \item{s}{MCMC chain for cell-specific technical normalisation parameters 
#' \eqn{s_j} (matrix with \code{n} columns, 
#' all elements must be positive numbers)}
#' \item{nu}{MCMC chain for cell-specific random effects \eqn{\nu_j}
#' (matrix with \code{n} columns, all elements must be positive numbers)}
#' \item{theta}{MCMC chain for technical over-dispersion parameter(s) 
#' \eqn{\theta} (matrix, all elements must be positive, 
#' each colum represents 1 batch)}
#' \item{\code{beta}}{MCMC chain for regression parameters
#' (matrix with \code{k} columns where k represent the number of chosen basis functions) }
#' \item{\code{sigma2}}{MCMC chain for the variance used during regression
#' (matrix with one column, sigma2 represents a gloabl parameter)}
#' \item{\code{eta}}{MCMC chain for the degrees of freedom used for regression
#' (matrix with one column)}
#' \item{\code{lambda}}{MCMC chain for the gene-specific random effect on the regression variance
#' (matrix with \code{q} columns)}
#' \item{\code{epsilon}}{MCMC chain for the gene specific regression residue (mean corrected vraribility)
#' (matrix with \code{q} columns)}
#' }
#' 
#' @return An object of class \code{\linkS4class{BASiCS_Chain}}.
#'
#' @examples
#'
#' Data <- makeExampleBASiCS_Data()
#' Chain <- BASiCS_MCMC(Data, N = 50, Thin = 5, Burn = 5)
#'
#' ChainMu <- displayChainBASiCS(Chain, 'mu')
#' ChainDelta <- displayChainBASiCS(Chain, 'delta')
#' ChainPhi <- displayChainBASiCS(Chain, 'phi')
#' ChainS <- displayChainBASiCS(Chain, 's')
#' ChainNu <- displayChainBASiCS(Chain, 'nu')
#' ChainTheta <- displayChainBASiCS(Chain, 'theta')
#' ChainBeta <- displayChainBASiCS(Chain, 'beta')
#' ChainSigma2 <- displayChainBASiCS(Chain, 'sigma2')
#' ChainEta <- displayChainBASiCS(Chain, 'eta')
#' ChainLambda <- displayChainBASiCS(Chain, 'lambda')
#' ChainEpsilon <- displayChainBASiCS(Chain, 'epsilon')
#'
#' ChainNew <- newBASiCS_Chain(parameters = list(mu = ChainMu, delta = ChainDelta,
#'                             phi = ChainPhi, s = ChainS, 
#'                             nu = ChainNu, theta = ChainTheta,
#'                             beta = ChainBeta, sigma2 = ChainSigma2,
#'                             eta = ChainEta, lambda = ChainLambda,
#'                             epsilon = ChainEpsilon))
#'
#' @seealso \code{\link[BASiCS]{BASiCS_Chain}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}, Nils Eling
#'
#' @references 
#' 
#' Vallejos, Marioni and Richardson (2015). PLoS Computational Biology. 
newBASiCS_Chain <- function(parameters) 
{
  Chain <- new("BASiCS_Chain", parameters = parameters)
  return(Chain)
}
