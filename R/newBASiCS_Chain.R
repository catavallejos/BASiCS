#' @title Creates a BASiCS_Chain object from pre-computed MCMC chains
#'
#' @description \code{BASiCS_Chain} creates a \code{\linkS4class{BASiCS_Chain}} 
#' object from pre-computed MCMC chains.
#'
#' @param parameters List of matrices containing MCMC chains for each model 
#' parameter. 
#' \describe{
##' \item{mu}{MCMC chain for gene-specific mean expression parameters 
#' \eqn{\mu_i}, biological genes only 
#' (matrix with \code{q.bio} columns, all elements must be positive numbers)}
#' \item{delta}{MCMC chain for gene-specific biological over-dispersion 
#' parameters \eqn{\delta_i}, biological genes only 
#' (matrix with \code{q.bio} columns, all elements must be positive numbers)}
#' \item{phi}{MCMC chain for cell-specific mRNA content normalisation parameters 
#' \eqn{\phi_j} (matrix with \code{n} columns, all elements must be positive 
#' numbers and the sum of its elements must be equal to \code{n})}
#' \item{s}{MCMC chain for cell-specific technical normalisation parameters 
#' \eqn{s_j} (matrix with \code{n} columns, 
#' all elements must be positive numbers)}
#' \item{nu}{MCMC chain for cell-specific random effects \eqn{\nu_j}
#' (matrix with \code{n} columns, all elements must be positive numbers)}
#' \item{theta}{MCMC chain for technical over-dispersion parameter(s) 
#' \eqn{\theta} (matrix, all elements must be positive, 
#' each colum represents 1 batch)}}
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
#'
#' ChainNew <- newBASiCS_Chain(parameters = list(mu = ChainMu, delta = ChainDelta,
#'                             phi = ChainPhi, s = ChainS, 
#'                             nu = ChainNu, theta = ChainTheta))
#'
#' @seealso \code{\link[BASiCS]{BASiCS_Chain}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#'
#' @references 
#' 
#' Vallejos, Marioni and Richardson (2015). PLoS Computational Biology. 
newBASiCS_Chain <- function(parameters) 
{
  Chain <- new("BASiCS_Chain", parameters = parameters)
  return(Chain)
}
