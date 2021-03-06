#' @title Generate a draw from the posterior of BASiCS using the
#' generative model.
#'
#' @description \code{BASiCS_Draw} creates a simulated dataset from the
#' posterior of a fitted model implemented in BASiCS.
#'
#' @param Chain An object of class \code{\linkS4class{BASiCS_Chain}}.
#' @param BatchInfo Vector of batch information from the SingleCellExperiment
#' object used as input to BASiCS_MCMC.
#' @param N The integer index for the draw to be used to sample from the 
#' posterior predictive distribution. If not supplied, a random value is chosen.
#' 
#' @return An object of class \code{\linkS4class{SingleCellExperiment}},
#' including synthetic data generated by the model implemented in BASiCS.
#'
#' @examples
#' data(ChainSC)
#' BASiCS_Draw(ChainSC)
#'
#' data(ChainSC)
#' BASiCS_Draw(ChainSC)
#'
#' @author Alan O'Callaghan
#'
#' @references
#'
#' Vallejos, Marioni and Richardson (2015). PLoS Computational Biology.
#'
#' @export
BASiCS_Draw <- function(
    Chain,
    BatchInfo = gsub(
      ".*_Batch([0-9a-zA-Z])", 
      "\\1",
      colnames(Chain@parameters[["nu"]])
    ),
    N = sample(nrow(Chain@parameters[["nu"]]), 1)
  ) {
  
  BASiCS_Sim(
    Mu = Chain@parameters[["mu"]][N, ],
    Delta = Chain@parameters[["delta"]][N, ],
    Phi = Chain@parameters[["phi"]][N, ],
    S = Chain@parameters[["s"]][N, ],
    BatchInfo = BatchInfo,
    Theta = Chain@parameters[["theta"]][N, ]
  )
}
