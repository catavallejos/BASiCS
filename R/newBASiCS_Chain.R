
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
#' \item{\code{beta}}{Only used for regression model. MCMC chain for regression 
#' coefficients (matrix with \code{k} columns, where \code{k} represent the 
#' number of chosen basis functions + 2) }
#' \item{\code{sigma2}}{Only used for regression model. MCMC chain for the 
#' residual variance (matrix with one column, sigma2 represents 
#' a global parameter)}
#' \item{\code{epsilon}}{Only used for regression model. MCMC chain for the 
#' gene specific residual over-dispersion parameter (mean corrected vraribility)
#' (matrix with \code{q} columns)}
#' }
#' 
#' @return An object of class \code{\linkS4class{BASiCS_Chain}}.
#'
#' @examples
#'
#' Data <- makeExampleBASiCS_Data()
#' 
#' # No regression model
#' Chain <- BASiCS_MCMC(Data, N = 50, Thin = 5, Burn = 5, Regression = FALSE)
#'
#' ChainMu <- displayChainBASiCS(Chain, 'mu')
#' ChainDelta <- displayChainBASiCS(Chain, 'delta')
#' ChainPhi <- displayChainBASiCS(Chain, 'phi')
#' ChainS <- displayChainBASiCS(Chain, 's')
#' ChainNu <- displayChainBASiCS(Chain, 'nu')
#' ChainTheta <- displayChainBASiCS(Chain, 'theta')
#'
#' ChainNew <- newBASiCS_Chain(parameters = list(mu = ChainMu, 
#'                                               delta = ChainDelta,
#'                                               phi = ChainPhi, 
#'                                               s = ChainS, 
#'                                               nu = ChainNu, 
#'                                               theta = ChainTheta))
#'                             
#'
#' # No regression model
#' Chain <- BASiCS_MCMC(Data, N = 50, Thin = 5, Burn = 5, Regression = TRUE)
#' 
#' ChainMu <- displayChainBASiCS(Chain, 'mu')
#' ChainDelta <- displayChainBASiCS(Chain, 'delta')
#' ChainPhi <- displayChainBASiCS(Chain, 'phi')
#' ChainS <- displayChainBASiCS(Chain, 's')
#' ChainNu <- displayChainBASiCS(Chain, 'nu')
#' ChainTheta <- displayChainBASiCS(Chain, 'theta')
#' ChainBeta <- displayChainBASiCS(Chain, 'beta')
#' ChainSigma2 <- displayChainBASiCS(Chain, 'sigma2')
#' ChainEpsilon <- displayChainBASiCS(Chain, 'epsilon')
#' 
#' ChainNew <- newBASiCS_Chain(parameters = list(mu = ChainMu, 
#'                                               delta = ChainDelta,
#'                                               phi = ChainPhi, 
#'                                               s = ChainS, 
#'                                               nu = ChainNu, 
#'                                               theta = ChainTheta,
#'                                               beta = ChainBeta, 
#'                                               sigma2 = ChainSigma2,
#'                                               epsilon = ChainEpsilon))
#'
#' @seealso \code{\link[BASiCS]{BASiCS_Chain}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#'
#' @export
newBASiCS_Chain <- function(parameters) {
  SampleParameters <- intersect(
    c("mu", "delta", "epsilon", "nu", "s", "theta", "phi", "beta", "sigma2", "lambda"),
    names(parameters)
  )
  parameters[SampleParameters] <- lapply(
    SampleParameters,
    function(Parameter) {
      x <- parameters[[Parameter]]
      attr(x, "ESS") <- ess(x)
      x
    }
  )
  Chain <- new("BASiCS_Chain", parameters = parameters)
  return(Chain)
}
