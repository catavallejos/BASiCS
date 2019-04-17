#' @name BASiCS_TestDE
#'
#' @title Remove global mean expression offset
#'
#' @description Remove global offset in mean expression between two 
#' \code{BASiCS_Chain} objects.
#'
#' @examples
#'
#' # Loading two 'BASiCS_Chain' objects (obtained using 'BASiCS_MCMC')
#' data(ChainSC)
#' data(ChainRNA)
#' 
#' BASiCS_CorrectOffset(ChainSC, ChainRNA, "a", "b", Plot = FALSE)
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#' @author Alan O'Callaghan \email{a.b.o'callaghan@sms.ed.ac.uk}
#' 
#' @export
BASiCS_CorrectOffset <- function(Chain1, 
                                 Chain2, 
                                 GroupLabel1 = "Group1", 
                                 GroupLabel2 = "Group2", 
                                 Plot = TRUE,
                                 ...) {
  
  n1 <- ncol(Chain1@parameters$nu)
  n2 <- ncol(Chain2@parameters$nu)
  n <- n1 + n2

  # Calculating iteration-specific offset
  OffsetChain <- matrixStats::rowSums2(Chain1@parameters$mu) /
                  matrixStats::rowSums2(Chain2@parameters$mu)
  # Offset point estimate
  OffsetEst <- median(OffsetChain)

  # Offset correction
  Chain1_offset <- Chain1
  Chain1_offset@parameters$mu <- Chain1@parameters$mu / OffsetEst
  Chain2_offset <- Chain2  # Chain2 requires no change
  Mu1 <- matrixStats::colMedians(Chain1_offset@parameters$mu)
  Mu2 <- matrixStats::colMedians(Chain2_offset@parameters$mu)
  Delta1 <- matrixStats::colMedians(Chain1_offset@parameters$delta)
  Delta2 <- matrixStats::colMedians(Chain2_offset@parameters$delta)

  Mu1_old <- matrixStats::colMedians(Chain1@parameters$mu)
  MuBase_old <- (Mu1_old * n1 + Mu2 * n2) / n
  ChainTau_old <- log2(Chain1@parameters$mu / Chain2@parameters$mu)
  MedianTau_old <- matrixStats::colMedians(ChainTau_old)

  # Offset corrected LFC estimates
  MuBase <- (Mu1 * n1 + Mu2 * n2)/n
  ChainTau <- log2(Chain1_offset@parameters$mu / Chain2_offset@parameters$mu)
  MedianTau <- matrixStats::colMedians(ChainTau)

  Corrected <- new("BASiCS_OffsetCorrected", 
    GroupLabel1 = GroupLabel1,
    GroupLabel2 = GroupLabel2,
    OffsetChain = OffsetChain,
    OffsetEst = OffsetEst,
    Chain1_offset = Chain1_offset,
    Chain2_offset = Chain2_offset,
    Mu1 = Mu1,
    Mu1_old = Mu1_old,
    Mu2 = Mu2,
    MuBase = MuBase,
    MuBase_old = MuBase_old,
    ChainTau = ChainTau,
    MedianTau = MedianTau,
    MedianTau_old = MedianTau_old,
    Delta1 = Delta1,
    Delta2 = Delta2
  )


  if (Plot) {
    BASiCS_PlotOffset(Corrected, ...)
  } else {
    message("-------------------------------------------------------------\n",
            "Offset estimate: ", round(OffsetEst, 4), "\n",
            "(ratio ", GroupLabel1, " vs ", GroupLabel2, ").\n",
            "To visualise its effect, please use 'PlotOffset = TRUE'.\n",
            "-------------------------------------------------------------\n")
  }
  Corrected
}
