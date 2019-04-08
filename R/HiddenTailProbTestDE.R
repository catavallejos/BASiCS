HiddenTailProbTestDE <- function(ChainLFC, Epsilon) {
  if (Epsilon > 0) {
    Prob <- matrixStats::colMeans2(ifelse(abs(ChainLFC) > Epsilon, 1, 0))
  }
  else {
    Prob_aux <- matrixStats::colMeans2(ifelse(ChainLFC > 0, 1, 0))
    Prob <- 2 * pmax(Prob_aux, 1 - Prob_aux) - 1
  }
  return(Prob)
}