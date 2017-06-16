#' @title Creates a BASiCS_Chain object from pre-computed MCMC chains
#'
#' @description \code{BASiCS_Chain} creates a \code{\link[BASiCS]{BASiCS_Chain-class}} object from pre-computed MCMC chains.
#'
#' @param mu MCMC chain for gene-specific expression levels \eqn{\mu[i]}, defined as true input molecules in case of technical genes
#' (matrix with \code{q} columns, technical genes located at the end of the matrix, all elements must be positive numbers)
#' @param delta MCMC chain for gene-specific biological cell-to-cell heterogeneity hyper-parameters \eqn{\delta[i]}, biological genes only
#' (matrix with \code{q.bio} columns, all elements must be positive numbers)
#' @param phi MCMC chain for cell-specific mRNA content normalising constants \eqn{\phi[j]}
#' (matrix with \code{n} columns, all elements must be positive numbers and the sum of its elements must be equal to \code{n})
#' @param s MCMC chain for cell-specific capture efficiency (or amplification biases if not using UMI based counts) normalising constants \eqn{s[j]}
#' (matrix with \code{n} columns, all elements must be positive numbers)
#' @param nu MCMC chain for cell-specific random effects \eqn{\nu[j]}
#' (matrix with \code{n} columns, all elements must be positive numbers)
#' @param theta MCMC chain for technical variability hyper-parameter \eqn{\theta} (vector, all elements must be positive)
#'
#' @return An object of class \code{\link[BASiCS]{BASiCS_Chain-class}}.
#'
#' @examples
#'
#' # Data = makeExampleBASiCS_Data()
#' # MCMC_Output <- BASiCS_MCMC(Data, N = 50, Thin = 5, Burn = 5,
#' #                StoreChains = TRUE, StoreDir = getwd(), RunName = "Test")
#'
#' # ChainMu = as.matrix(read.table("chain_mu_Test.txt"))
#' # ChainDelta = as.matrix(read.table("chain_delta_Test.txt"))
#' # ChainPhi = as.matrix(read.table("chain_phi_Test.txt"))
#' # ChainS = as.matrix(read.table("chain_s_Test.txt"))
#' # ChainNu = as.matrix(read.table("chain_nu_Test.txt"))#
#' # ChainTheta = read.table("chain_theta_Test.txt")[,1]
#'
#' # MCMC_Output_Load <- newBASiCS_Chain(mu = ChainMu, delta = ChainDelta,
#' #   phi = ChainPhi, s = ChainS, nu = ChainNu, theta = ChainTheta)
#'
#' @seealso \code{\link[BASiCS]{BASiCS_Chain-class}}
#'
#' @author Catalina A. Vallejos \email{cvallejos@turing.ac.uk}
#'
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
newBASiCS_Chain <- function(mu, delta, phi, s, nu, theta)
{
  Chain <- new("BASiCS_Chain", mu = mu, delta = delta, phi = phi, s = s, nu = nu, theta = theta)
  show(Chain)
  return(Chain)
}