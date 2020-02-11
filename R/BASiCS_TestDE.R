#' @name BASiCS_TestDE
#'
#' @title Detection of genes with changes in expression
#'
#' @description Function to assess changes in expression between two groups
#' of cells (mean and over-dispersion)
#'
#' @param Chain1 an object of class \code{\linkS4class{BASiCS_Chain}}
#' containing parameter estimates for the first group of cells
#' @param Chain2 an object of class \code{\linkS4class{BASiCS_Chain}}
#' containing parameter estimates for the second group of cells
#' @param EpsilonM Minimum fold change tolerance threshold for detecting
#' changes in overall expression (must be a positive real number).
#' Default value: \code{EpsilonM = log2(1.5)} (i.e. 50\% increase).
#' @param EpsilonD Minimum fold change tolerance threshold for detecting
#' changes in biological over-dispersion (must be a positive real number).
#' Default value: \code{EpsilonM = log2(1.5)} (i.e. 50\% increase).
#' @param EpsilonR Minimum distance threshold for detecting
#' changes in residual over-dispersion (must be a positive real number).
#' Default value: \code{EpsilonR= log2(1.5)/log2(exp(1))} (i.e. 50\% increase).
#' @param ProbThresholdM Optional parameter. Probability threshold for detecting
#' changes in overall expression (must be a positive value, between 0 and 1).
#' If \code{EFDR_M = NULL}, the posterior probability threshold for the
#' differential mean expression test will be set to \code{ProbThresholdM}. If
#' a value for \code{EFDR_M} is provided, the posterior probability threshold
#' is chosen to achieve an EFDR equal to \code{EFDR_M} and \code{ProbThresholdM}
#' defines a minimum probability threshold for this calibration (this avoids low
#' values of \code{ProbThresholdM} to be chosen by the EFDR calibration.
#' Default value \code{ProbThresholdM = 2/3}, i.e. the probability of observing
#' a log2-FC above \code{EpsilonM} must be at least twice the probality of
#' observing the complementary event (log2-FC below \code{EpsilonM}).
#' @param ProbThresholdD Optional parameter. Probability threshold for detecting
#' changes in cell-to-cell biological over-dispersion (must be a positive value,
#' between 0 and 1). Same usage as \code{ProbThresholdM}, depending on the value
#' provided for \code{EFDR_D}. Default value \code{ProbThresholdD = 2/3}.
#' @param ProbThresholdR Optional parameter. Probability threshold for detecting
#' changes in residual over-dispersion (must be a positive value, between 0 and
#' 1). Same usage as \code{ProbThresholdM}, depending on the value provided for
#' \code{EFDR_R}. Default value \code{ProbThresholdR = 2/3}.
#' @param OrderVariable Ordering variable for output.
#' Possible values: \code{'GeneIndex'} (default), \code{'GeneName'} and
#' \code{'Mu'} (mean expression).
#' @param GroupLabel1 Label assigned to reference group.
#' Default: \code{GroupLabel1 = 'Group1'}
#' @param GroupLabel2 Label assigned to reference group.
#' Default: \code{GroupLabel2 = 'Group2'}
#' @param Plot If \code{Plot = TRUE}, MA and volcano plots are generated.
#' @param PlotOffset If \code{Plot = TRUE}, the offset effect is visualised.
#' @param Offset Optional argument to remove a fix offset effect (if not
#' previously removed from the MCMC chains). Default: \code{Offset = TRUE}.
#' @param EFDR_M Target for expected false discovery rate related to
#' the comparison of means. If \code{EFDR_M = NULL}, EFDR calibration is not
#' performed and the posterior probability threshold is set equal to
#' \code{ProbThresholdM}. Default \code{EFDR_M = 0.05}.
#' @param EFDR_D Target for expected false discovery rate related to
#' the comparison of dispersions. If \code{EFDR_D = NULL}, EFDR calibration is
#' not performed and the posterior probability threshold is set equal to
#' \code{ProbThresholdD}.Default \code{EFDR_D = 0.05}.
#' @param EFDR_R Target for expected false discovery rate related to
#' the comparison of residual over-dispersions. If \code{EFDR_R = NULL}, EFDR
#' calibration is not performed and the posterior probability threshold is set
#' equal to \code{ProbThresholdR}.Default \code{EFDR_D = 0.05}.
#' @param GenesSelect Optional argument to provide a user-defined list
#' of genes to be considered for the comparison.
#' Default: \code{GenesSelect = rep(TRUE, nGene)}
#' When used, this argument must be a vector
#' of \code{TRUE} (include gene) / \code{FALSE} (exclude gene) indicator,
#' with the same length as the number of intrinsic genes and following the same
#' order as how genes are displayed in the table of counts.
#' This argument is necessary in order to have a meaningful EFDR calibration
#' when the user decides to exclude some genes from the comparison.
#' @param min.mean Minimum mean expression threshold required for inclusion in
#' offset calculation. Similar to `min.mean` in `scran::computeSumFactors`. This
#' parameter is only relevant with `Offset = TRUE`.
#' @param CheckESS Should the effective sample size of the chains be tested to
#' be of a suitable magnitude in order to be included in tests for
#' differential expression? Default is FALSE. If TRUE, genes with poor
#' mixing will be excluded from the tests for differential expression.
#' @param ESSThreshold If CheckESS is TRUE, this argument specifies
#' the minimum effective sample size for a gene to be included in the tests for
#' differential expression. Default is 100.
#' @param ... Optional parameters.
#'
#' @return \code{BASiCS_TestDE} returns a list of 4 elements:
#' \describe{
#' \item{\code{TableMean}}{A \code{\link[base]{data.frame}} containing the
#'       results of the differential mean test
#'    \describe{
#'    \item{\code{GeneName}}{Gene name}
#'    \item{\code{MeanOverall}}{For each gene, the estimated mean expression
#'          parameter \eqn{\mu_i} is averaged across both groups of cells
#'          (weighted by sample size).}
#'    \item{\code{Mean1}}{Estimated mean expression parameter \eqn{\mu_i} for
#'          each biological gene in the first group of cells.}
#'    \item{\code{Mean2}}{Estimated mean expression parameter \eqn{\mu_i} for
#'          each biological gene in the second group of cells.}
#'    \item{\code{MeanFC}}{Fold change in mean expression parameters between
#'          the first and second groups of cells.}
#'    \item{\code{MeanLog2FC}}{Log2-transformed fold change in mean expression
#'          between the first and second groups of cells.}
#'    \item{\code{ProbDiffMean}}{Posterior probability for mean expression
#'          difference between the first and second groups of cells.}
#'    \item{\code{ResultDiffExp}}{Indicator if a gene has a higher mean
#'          expression in the first or second groups of cells.}
#'    }}
#' \item{\code{TableDisp}}{A \code{\link[base]{data.frame}} containing the
#'       results of the differential dispersion test (excludes genes for which
#'       the mean does not changes).
#'    \describe{
#'    \item{\code{GeneName}}{Gene name}
#'    \item{\code{MeanOverall}}{For each gene, the estimated mean expression
#'          parameter \eqn{\mu_i} is averaged across both groups of cells
#'          (weighted by sample size).}
#'    \item{\code{DispOverall}}{For each gene, the estimated over-dispersion
#'          parameter \eqn{\delta_i} is averaged across both groups of cells
#'          (weighted by sample size).}
#'    \item{\code{Disp1}}{Estimated over-dispersion parameter \eqn{\delta_i}
#'          for each biological gene in the first group of cells.}
#'    \item{\code{Disp2}}{Estimated over-dispersion parameter \eqn{\delta_i}
#'          for each biological gene in the second group of cells.}
#'    \item{\code{DispFC}}{Fold change in over-dispersion parameters between
#'          the between the first and second groups of cells.}
#'    \item{\code{DispLog2FC}}{Log-transformed fold change in over-dispersion
#'          between the first and second groups of cells.}
#'    \item{\code{ProbDiffDisp}}{Posterior probability for over-dispersion
#'          difference between the first and second groups of cells.}
#'    \item{\code{ResultDiffDisp}}{Indicator if a gene has a higher
#'          over-dispersion in the first or second groups of cells. Genes
#'          labelled with "ExcludedFromTest" were detected as showing
#'          differential mean expression.}
#'    }}
#' \item{\code{TableResDisp}}{A \code{\link[base]{data.frame}} containing the
#'       results of the differential residual over-dispersion test.
#'    \describe{
#'    \item{\code{GeneName}}{Gene name}
#'    \item{\code{MeanOverall}}{For each gene, the estimated mean expression
#'          parameter \eqn{\mu_i} is averaged across both groups of cells
#'          (weighted by sample size).}
#'    \item{\code{ResDispOverall}}{For each gene, the estimated residual
#'          over-dispersion parameter \eqn{\delta_i} is averaged across both
#'          groups of cells (weighted by sample size).}
#'    \item{\code{ResDisp1}}{Estimated residual over-dispersion parameter
#'          \eqn{\epsilon_i} for each biological gene in the first group of
#'          cells.}
#'    \item{\code{ResDisp2}}{Estimated residual over-dispersion parameter
#'          \eqn{\epsilon_i} for each biological gene in the second group of
#'          cells.}
#'    \item{\code{ResDispDistance}}{Difference in residual over-dispersion
#'          between the first and second groups of cells.}
#'    \item{\code{ProbDiffResDisp}}{Posterior probability for residual
#'          over-dispersion difference between the first and second groups of
#'          cells.}
#'    \item{\code{ResultDiffResDisp}}{Indicator if a gene has a higher
#'          residual over-dispersion in the first or second groups of cells.
#'          Genes labelled with "ExcludedFromTest" were not expressed in
#'          at least 2 cells per condition.}
#'    }}
#'  \item{\code{DiffMeanSummary}}{A list containing the following information
#'        for the differential mean expression test:
#'    \describe{
#'   \item{\code{ProbThreshold}}{Posterior probability threshold.}
#'   \item{\code{EFDR}}{Expected false discovery rate for the given thresholds.}
#'   \item{\code{EFNR}}{Expected false negative rate for the given thresholds.}
#'   }}
#'  \item{\code{DiffDispSummary}}{A list containing the following
#'        information for the differential over-dispersion test:
#'    \describe{
#'   \item{\code{ProbThreshold}}{Posterior probability threshold.}
#'   \item{\code{EFDR}}{Expected false discovery rate for the given thresholds.}
#'   \item{\code{EFNR}}{Expected false negative rate for the given thresholds.}
#'   }}
#'  \item{\code{DiffResDispSummary}}{A list containing the following
#'        information for the differential residual over-dispersion test:
#'    \describe{
#'   \item{\code{ProbThreshold}}{Posterior probability threshold.}
#'   \item{\code{EFDR}}{Expected false discovery rate for the given thresholds.}
#'   \item{\code{EFNR}}{Expected false negative rate for the given thresholds.}
#'   }}
#'  \item{\code{Chain1_offset}}{an \code{\link[BASiCS]{BASiCS_Chain}} object:
#'        \code{Chain1} after offset removal.}
#'  \item{\code{Chain2_offset}}{an \code{\link[BASiCS]{BASiCS_Chain}} object:
#'        \code{Chain2} after offset removal (this is only provided for
#'        completeness; \code{Chain2} is not affected by the offset).}
#'  \item{\code{OffsetChain}}{MCMC chain calculated for the offset effect.}
#'  \item{\code{Offset}}{Estimated offset (posterior median of
#'        \code{OffsetChain}). Default value set equal to 1 when offset
#'        correction is not performed.}
#' }
#'
#' @examples
#'
#' # Loading two 'BASiCS_Chain' objects (obtained using 'BASiCS_MCMC')
#' data(ChainSC)
#' data(ChainRNA)
#'
#' Test <- BASiCS_TestDE(
#'   Chain1 = ChainSC, Chain2 = ChainRNA,
#'   GroupLabel1 = "SC", GroupLabel2 = "P&S",
#'   EpsilonM = log2(1.5), EpsilonD = log2(1.5),
#'   OffSet = TRUE
#' )
#'
#' # Results for the differential mean test
#' head(Test$TableMean)
#'
#' # Results for the differential over-dispersion test
#' # This only includes genes marked as 'NoDiff' in Test$TableMean
#' head(Test$TableDisp)
#'
#' # For testing differences in residual over-dispersion, two chains obtained
#' # via 'BASiCS_MCMC(Data, N, Thin, Burn, Regression=TRUE)' need to be provided
#' data(ChainSCReg)
#' data(ChainRNAReg)
#'
#' Test <- BASiCS_TestDE(
#'   Chain1 = ChainSCReg, Chain2 = ChainRNAReg,
#'   GroupLabel1 = "SC", GroupLabel2 = "P&S",
#'   EpsilonM = log2(1.5), EpsilonD = log2(1.5),
#'   EpsilonR = log2(1.5) / log2(exp(1)),
#'   OffSet = TRUE
#' )
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#'
#' @rdname BASiCS_TestDE
#' @export
BASiCS_TestDE <- function(Chain1,
                          Chain2,
                          EpsilonM = log2(1.5),
                          EpsilonD = log2(1.5),
                          EpsilonR = log2(1.5) / log2(exp(1)),
                          ProbThresholdM = 2 / 3,
                          ProbThresholdD = 2 / 3,
                          ProbThresholdR = 2 / 3,
                          OrderVariable = "GeneIndex",
                          GroupLabel1 = "Group1",
                          GroupLabel2 = "Group2",
                          Plot = TRUE,
                          PlotOffset = TRUE,
                          Offset = TRUE,
                          EFDR_M = 0.05,
                          EFDR_D = 0.05,
                          EFDR_R = 0.05,
                          GenesSelect = rep(TRUE, ncol(Chain1@parameters[["mu"]])),
                          min.mean = 1,
                          CheckESS = FALSE,
                          ESSThreshold = 100,
                          ...) {
  HiddenHeaderTest_DE(
    Chain1 = Chain1,
    Chain2 = Chain2,
    EpsilonM = EpsilonM,
    EpsilonD = EpsilonD,
    EpsilonR = EpsilonR,
    EFDR_M = EFDR_M,
    EFDR_D = EFDR_D,
    EFDR_R = EFDR_R,
    ProbThresholdM = ProbThresholdM,
    ProbThresholdD = ProbThresholdD,
    ProbThresholdR = ProbThresholdR,
    OrderVariable = OrderVariable,
    GroupLabel1 = GroupLabel1,
    GroupLabel2 = GroupLabel2,
    GenesSelect = GenesSelect,
    Plot = Plot,
    PlotOffset = PlotOffset,
    Offset = Offset,
    CheckESS = CheckESS
  )

  GeneName <- colnames(Chain1@parameters$mu)
  GeneIndex <- seq_len(length(GeneName))

  # If all genes are to be included
  if (is.null(GenesSelect)) GenesSelect <- rep(TRUE, times = length(GeneName))

  message(
    "-------------------------------------------------------------\n",
    "Log-fold change thresholds are now set in a log2 scale. \n",
    "Original BASiCS release used a natural logarithm scale."
  )


  if (xor(is.null(Chain1@parameters[["epsilon"]]), is.null(Chain2@parameters[["epsilon"]]))) {
    stop("Both chains should be run with the same setting for Regression.")
  }

  n1 <- ncol(Chain1@parameters$nu)
  n2 <- ncol(Chain2@parameters$nu)
  n <- n1 + n2
  # With offset correction
  if (Offset) {
    # Offset correction
    A <- BASiCS_CorrectOffset(Chain1, Chain2, min.mean = min.mean)
    OffsetEst <- A$Offset
    OffsetChain <- A$OffsetChain
    Chain1_offset <- A$Chain
    Chain2_offset <- Chain2 # Chain2 requires no change

    # Post-offset correction gene-specific estimates
    Mu1 <- matrixStats::colMedians(Chain1_offset@parameters$mu)
    Mu2 <- matrixStats::colMedians(Chain2_offset@parameters$mu)
    Delta1 <- matrixStats::colMedians(Chain1_offset@parameters$delta)
    Delta2 <- matrixStats::colMedians(Chain2_offset@parameters$delta)

    # Pre-offset correction LFC estimates
    Mu1_old <- matrixStats::colMedians(Chain1@parameters$mu)
    MuBase_old <- (Mu1_old * n1 + Mu2 * n2) / n
    ChainTau_old <- log2(Chain1@parameters$mu / Chain2@parameters$mu)
    MedianTau_old <- matrixStats::colMedians(ChainTau_old)

    # Offset corrected LFC estimates
    MuBase <- (Mu1 * n1 + Mu2 * n2) / n
    ChainTau <- log2(Chain1_offset@parameters$mu / Chain2_offset@parameters$mu)
    MedianTau <- matrixStats::colMedians(ChainTau)

    if (!PlotOffset) {
      message(
        "-------------------------------------------------------------\n",
        "Offset estimate: ", round(OffsetEst, 4), "\n",
        "(ratio ", GroupLabel1, " vs ", GroupLabel2, ").\n",
        "To visualise its effect, please use 'PlotOffset = TRUE'.\n",
        "-------------------------------------------------------------\n"
      )
    } else {
      message(
        "-------------------------------------------------------------\n",
        "Offset estimate: ", round(OffsetEst, 4), "\n",
        "(ratio ", GroupLabel1, " vs ", GroupLabel2, ").\n",
        "-------------------------------------------------------------\n"
      )
    }

    if (PlotOffset) {
      message(
        "Plots to follow: \n",
        "1. Posterior uncertainty for offset estimate \n",
        "2. Mean expression estimates before/after offset correction \n",
        "3. MA plot before/after offset correction \n"
      )

      par(ask = TRUE)
      # Offset uncertainty
      graphics::boxplot(OffsetChain,
        frame = FALSE,
        main = "Offset MCMC chain", ylab = "Offset estimate"
      )
      # Mean expression parameters before/after offset correction
      par(mfrow = c(1, 2))
      graphics::boxplot(cbind(Mu1_old, Mu2),
        frame = FALSE, main = "Before correction",
        names = c(GroupLabel1, GroupLabel2),
        ylab = "Mean expression", log = "y"
      )
      graphics::boxplot(cbind(Mu1, Mu2),
        frame = FALSE, main = "After correction",
        names = c(GroupLabel1, GroupLabel2),
        ylab = "Mean expression", log = "y"
      )
      # MA plot pre/after offset
      par(mfrow = c(1, 2))
      graphics::smoothScatter(log2(MuBase_old), MedianTau_old,
        bty = "n",
        xlab = "Mean expresssion (log2)",
        ylab = paste(
          "Log2 fold change", GroupLabel1,
          "vs", GroupLabel2
        ),
        main = "Before correction"
      )
      abline(h = 0, lty = 2)
      abline(h = log2(OffsetEst), lty = 1, col = "red")
      legend("topright", "log2offset", lty = 1, col = "red")
      graphics::smoothScatter(log2(MuBase), MedianTau,
        bty = "n",
        xlab = "Mean expresssion (log2)",
        ylab = paste(
          "Log2 fold change", GroupLabel1,
          "vs", GroupLabel2
        ),
        main = "After correction"
      )
      abline(h = 0, lty = 2)
      par(ask = FALSE)
    }
  } else {
    message(
      "-------------------------------------------------------------\n",
      "It is recomended to perform a global offset correction \n",
      "to remove global changes between the two groups of cells \n",
      "Default offset value set equal to 1.\n",
      "To perform offset correction, please set 'Offset = TRUE'. \n",
      "-------------------------------------------------------------\n"
    )

    # Summary1 <- Summary(Chain1)
    # Summary2 <- Summary(Chain2)
    Mu1 <- matrixStats::colMedians(Chain1@parameters$mu)
    Mu2 <- matrixStats::colMedians(Chain2@parameters$mu)
    Delta1 <- matrixStats::colMedians(Chain1@parameters$delta)
    Delta2 <- matrixStats::colMedians(Chain2@parameters$delta)
    MuBase <- (Mu1 * n1 + Mu2 * n2) / n
    ChainTau <- log2(Chain1@parameters$mu / Chain2@parameters$mu)
    MedianTau <- matrixStats::colMedians(ChainTau)

    # Default values when no offset correction is applied
    OffsetEst <- 1
    OffsetChain <- NULL
    Chain1_offset <- NULL
    Chain2_offset <- NULL
  }

  Search <- is.null(ProbThresholdM)

  if (CheckESS) {
    GoodESS <- ess(mcmc(Chain1@parameters[["mu"]])) > ESSThreshold &
      ess(mcmc(Chain2@parameters[["mu"]])) > ESSThreshold
  } else {
    GoodESS <- rep(TRUE, length(GenesSelect))
  }

  # Changes in mean expression
  # Calculating posterior probabilities
  ProbM <- .TailProb(Chain = abs(ChainTau), Threshold = EpsilonM)
  AuxMean <- .ThresholdSearch(
    Probs = ProbM[GenesSelect],
    ProbThreshold = ProbThresholdM,
    EFDR = EFDR_M,
    Task = "Differential mean",
    Suffix = "M"
  )
  OptThresholdM <- AuxMean$OptThreshold
  # Test results
  ResultDiffMean <- .TestResults(
    Probs = ProbM,
    Threshold = OptThresholdM[1],
    Estimate = MedianTau,
    Label1 = GroupLabel1,
    Label2 = GroupLabel2,
    GenesSelect = GenesSelect,
    GoodESS = GoodESS
  )

  # Output table
  TableMean <- cbind.data.frame(
    GeneName = GeneName,
    MeanOverall = as.numeric(MuBase),
    Mean1 = Mu1,
    Mean2 = Mu2,
    MeanFC = as.numeric(2^(MedianTau)),
    MeanLog2FC = as.numeric(MedianTau),
    ProbDiffMean = as.numeric(ProbM),
    ResultDiffMean = ResultDiffMean,
    stringsAsFactors = FALSE
  )
  # Rounding to 3 decimal points
  TableMean[, 2:7] <- round(TableMean[, 2:7], 3)

  # Genes with no change in mean expression
  NotDE <- ResultDiffMean == "NoDiff"

  # Changes in over dispersion
  ChainOmega <- log2(Chain1@parameters$delta / Chain2@parameters$delta)
  MedianOmega <- matrixStats::colMedians(ChainOmega)
  DeltaBase <- (Delta1 * n1 + Delta2 * n2) / n

  # Genes to calibrate EFDR
  if (!is.null(GenesSelect)) {
    DeltaSelect <- NotDE & GenesSelect
  } else {
    DeltaSelect <- NotDE
  }

  if (CheckESS) {
    GoodESS <- ess(mcmc(Chain1@parameters[["delta"]])) > ESSThreshold &
      ess(mcmc(Chain2@parameters[["delta"]])) > ESSThreshold
  }
  DeltaSelect <- DeltaSelect & GoodESS


  ProbD <- .TailProb(Chain = abs(ChainOmega), Threshold = EpsilonD)
  AuxDisp <- .ThresholdSearch(
    Probs = ProbD[DeltaSelect],
    ProbThreshold = ProbThresholdD,
    EFDR = EFDR_D,
    Task = "Differential dispersion",
    Suffix = "D"
  )
  OptThresholdD <- AuxDisp$OptThreshold

  # Test results
  ResultDiffDisp <- .TestResults(
    Probs = ProbD,
    Threshold = OptThresholdD[1],
    Estimate = MedianOmega,
    Label1 = GroupLabel1,
    Label2 = GroupLabel2,
    GenesSelect = GenesSelect,
    GoodESS = GoodESS,
    Excluded = !NotDE
  )

  # Output table
  TableDisp <- cbind.data.frame(
    GeneName = GeneName,
    MeanOverall = as.numeric(MuBase),
    DispOverall = as.numeric(DeltaBase),
    Disp1 = Delta1,
    Disp2 = Delta2,
    DispFC = as.numeric(2^(MedianOmega)),
    DispLog2FC = as.numeric(MedianOmega),
    ProbDiffDisp = as.numeric(ProbD),
    ResultDiffDisp = ResultDiffDisp,
    stringsAsFactors = FALSE
  )
  # Rounding to 3 decimal points
  TableDisp[, 2:8] <- round(TableDisp[, 2:8], 3)

  # Ordering the output tables
  if (OrderVariable == "GeneIndex") {
    orderVar <- order(GeneIndex, decreasing = FALSE)
  }
  if (OrderVariable == "GeneName") {
    orderVar <- order(GeneName, decreasing = TRUE)
  }
  if (OrderVariable == "Mu") {
    orderVar <- order(as.numeric(MuBase), decreasing = TRUE)
  }
  TableMean <- TableMean[orderVar, ]
  TableDisp <- TableDisp[orderVar, ]

  # Changes in residual over-dispersion - if regression approach was used
  if (!is.null(Chain1@parameters$epsilon)) {
    NotExcluded <- !(is.na(Chain1@parameters$epsilon[1, ]) |
      is.na(Chain2@parameters$epsilon[1, ]))

    ChainPsi <- Chain1@parameters$epsilon - Chain2@parameters$epsilon
    MedianPsi <- matrixStats::colMedians(ChainPsi)
    Epsilon1 <- matrixStats::colMedians(Chain1@parameters$epsilon)
    Epsilon2 <- matrixStats::colMedians(Chain2@parameters$epsilon)
    EpsilonBase <- (Epsilon1 * n1 + Epsilon2 * n2) / n

    # Genes to calibrate EFDR
    if (!is.null(GenesSelect)) {
      EpsSelect <- NotExcluded & GenesSelect
    } else {
      EpsSelect <- NotExcluded
    }
    if (CheckESS) {
      GoodESS <- ess(mcmc(Chain1@parameters[["epsilon"]])) > ESSThreshold &
        ess(mcmc(Chain2@parameters[["epsilon"]])) > ESSThreshold
    }
    EpsSelect <- EpsSelect & GoodESS

    select <- NotExcluded & GenesSelect
    ProbE <- .TailProb(Chain = abs(ChainPsi), Threshold = EpsilonR)
    AuxResDisp <- .ThresholdSearch(
      Probs = ProbE[EpsSelect],
      ProbThreshold = ProbThresholdR,
      EFDR = EFDR_R,
      Task = "Differential residual dispersion",
      Suffix = "R"
    )
    OptThresholdE <- AuxResDisp$OptThreshold

    # Test results
    ResultDiffResDisp <- .TestResults(
      Probs = ProbE,
      Threshold = OptThresholdE[1],
      Estimate = MedianPsi,
      Label1 = GroupLabel1,
      Label2 = GroupLabel2,
      GenesSelect = GenesSelect,
      GoodESS = GoodESS,
      Excluded = !NotExcluded
    )

    # Output table
    TableResDisp <- cbind.data.frame(
      GeneName = GeneName,
      MeanOverall = as.numeric(MuBase),
      ResDispOverall = as.numeric(EpsilonBase),
      ResDisp1 = Epsilon1,
      ResDisp2 = Epsilon2,
      ResDispDistance = as.numeric(MedianPsi),
      ProbDiffResDisp = as.numeric(ProbE),
      ResultDiffResDisp = ResultDiffResDisp,
      stringsAsFactors = FALSE
    )

    # Rounding to 3 decimal points
    TableResDisp[, 2:7] <- round(TableResDisp[, 2:7], 3)
    TableResDisp <- TableResDisp[orderVar, ]
  }



  if (!is.null(GenesSelect)) {
    message(
      "-------------------------------------------------------------\n",
      "The user excluded ", sum(!GenesSelect), " genes. \n",
      "These genes are marked as 'ExcludedByUser' \n",
      "and excluded from EFDR calibration. \n",
      "-------------------------------------------------------------\n"
    )
  }

  if (Plot) {
    if (Search) {
      message(
        "Plots to follow: \n",
        "1. EFDR/EFNR control plots \n",
        "2. MA plots \n",
        "3. Volcano plots \n"
      )
    } else {
      message(
        "Plots to follow: \n",
        "1. MA plots \n",
        "2. Volcano plots \n"
      )
    }

    par(ask = TRUE)

    if (Search) {
      if (!is.null(Chain1@parameters$epsilon)) {
        par(mfrow = c(1, 3))
      } else {
        par(mfrow = c(1, 2))
      }
      ProbThresholds <- seq(0.5, 0.9995, by = 0.00025)
      plot(ProbThresholds, AuxMean$EFDRgrid,
        type = "l", lty = 1, bty = "n",
        ylab = "Error rate", xlab = "Probability threshold",
        ylim = c(0, 1), main = "Differential mean"
      )
      lines(ProbThresholds, AuxMean$EFNRgrid, lty = 2)
      abline(h = EFDR_M, col = "blue", lwd = 2, lty = 1)
      abline(v = OptThresholdM[1], col = "red", lwd = 2, lty = 1)
      legend("top", c("EFDR", "EFNR", "Target EFDR"),
        lty = c(1, 2, 1),
        col = c("black", "black", "blue"), bty = "n"
      )
      plot(ProbThresholds, AuxDisp$EFDRgrid,
        type = "l", lty = 1, bty = "n",
        ylab = "Error rate", xlab = "Probability threshold",
        ylim = c(0, 1), main = "Differential dispersion"
      )
      lines(ProbThresholds, AuxDisp$EFNRgrid, lty = 2)
      abline(h = EFDR_D, col = "blue", lwd = 2, lty = 1)
      abline(v = OptThresholdD[1], col = "red", lwd = 2, lty = 1)
      legend("top", c("EFDR", "EFNR", "Target EFDR"),
        lty = c(1, 2, 1),
        col = c("black", "black", "blue"), bty = "n"
      )
      if (!is.null(Chain1@parameters$epsilon)) {
        plot(ProbThresholds, AuxResDisp$EFDRgrid,
          type = "l", lty = 1, bty = "n",
          ylab = "Error rate", xlab = "Probability threshold",
          ylim = c(0, 1), main = "Differential residual dispersion"
        )
        lines(ProbThresholds, AuxResDisp$EFNRgrid, lty = 2)
        abline(h = EFDR_R, col = "blue", lwd = 2, lty = 1)
        abline(v = OptThresholdE[1], col = "red", lwd = 2, lty = 1)
        legend("top", c("EFDR", "EFNR", "Target EFDR"),
          lty = c(1, 2, 1),
          col = c("black", "black", "blue"), bty = "n"
        )
      }
    }

    # MA plots
    if (!is.null(Chain1@parameters$epsilon)) {
      par(mfrow = c(1, 3))
    } else {
      par(mfrow = c(1, 2))
    }
    graphics::smoothScatter(
      log2(TableMean[["MeanOverall"]]),
      TableMean[["MeanLog2FC"]],
      bty = "n",
      xlab = "Mean expresssion (log2)",
      ylab = paste(
        "Log2 fold change",
        GroupLabel1, "vs",
        GroupLabel2
      ),
      main = "Differential mean"
    )
    indUse <- !(TableMean$ResultDiffMean %in% c("ExcludedByUser", "NoDiff"))
    points(
      log2(TableMean[indUse, "MeanOverall"]),
      TableMean[indUse, "MeanLog2FC"],
      pch = 16,
      col = "red"
    )
    
    abline(h = c(-EpsilonM, EpsilonM), lty = 2)
    
    graphics::smoothScatter(
      log2(TableDisp[["MeanOverall"]]),
      TableDisp[["DispLog2FC"]],
      bty = "n",
      xlab = "Mean expresssion (log2)",
      ylab = paste(
        "Log2 fold change",
        GroupLabel1, "vs",
        GroupLabel2
      ),
      main = "Differential dispersion"
    )

    with(
      TableDisp[!(TableDisp$ResultDiffDisp %in%
        c("ExcludedFromTesting", "ExcludedByUser", "NoDiff")), ],
      points(log2(MeanOverall), DispLog2FC, pch = 16, col = "red")
    )
    abline(h = c(-EpsilonD, EpsilonD), lty = 2)

    if (!is.null(Chain1@parameters$epsilon)) {
      with(
        TableResDisp[TableResDisp$ResultDiffResDisp != "ExcludedFromTesting", ],
        graphics::smoothScatter(log2(MeanOverall), ResDispDistance,
          bty = "n",
          xlab = "Mean expresssion (log2)",
          ylab = paste(
            "Difference",
            GroupLabel1, "vs",
            GroupLabel2
          ),
          main = "Differential residual dispersion"
        )
      )
      with(
        TableResDisp[!(TableResDisp$ResultDiffResDisp %in%
          c("ExcludedFromTesting", "ExcludedByUser", "NoDiff")), ],
        points(log2(MeanOverall), ResDispDistance, pch = 16, col = "red")
      )
      abline(h = c(-EpsilonR, EpsilonR), lty = 2)
    }

    # Volcano plots
    if (!is.null(Chain1@parameters$epsilon)) {
      par(mfrow = c(1, 3))
    } else {
      par(mfrow = c(1, 2))
    }
    with(
      TableMean,
      graphics::smoothScatter(MeanLog2FC, ProbDiffMean,
        bty = "n", ylim = c(0, 1),
        ylab = "Posterior probability",
        xlab = paste(
          "Log2 fold change",
          GroupLabel1, "vs",
          GroupLabel2
        ),
        main = "Differential mean test"
      )
    )
    with(
      TableMean[!(TableMean$ResultDiffMean %in%
        c("ExcludedByUser", "NoDiff")), ],
      points(MeanLog2FC, ProbDiffMean, pch = 16, col = "red")
    )
    abline(v = c(-EpsilonM, EpsilonM), lty = 2)
    with(
      TableDisp,
      graphics::smoothScatter(DispLog2FC, ProbDiffDisp,
        bty = "n", ylim = c(0, 1),
        ylab = "Posterior probability",
        xlab = paste(
          "Log2 fold change",
          GroupLabel1, "vs",
          GroupLabel2
        ),
        main = "Differential dispersion test"
      )
    )
    with(
      TableDisp[!(TableDisp$ResultDiffDisp %in%
        c("ExcludedFromTesting", "ExcludedByUser", "NoDiff")), ],
      points(DispLog2FC, ProbDiffDisp, pch = 16, col = "red")
    )
    abline(v = c(-EpsilonD, EpsilonD), lty = 2)

    if (!is.null(Chain1@parameters$epsilon)) {
      with(
        TableResDisp[TableResDisp$ResultDiffResDisp != "ExcludedFromTesting", ],
        graphics::smoothScatter(ResDispDistance, ProbDiffResDisp,
          bty = "n", ylim = c(0, 1),
          ylab = "Posterior probability",
          xlab = paste(
            "Difference",
            GroupLabel1, "vs",
            GroupLabel2
          ),
          main = "Differential residual dispersion test"
        )
      )
      with(
        TableResDisp[!(TableResDisp$ResultDiffResDisp %in%
          c("ExcludedFromTesting", "ExcludedByUser", "NoDiff")), ],
        points(ResDispDistance, ProbDiffResDisp, pch = 16, col = "red")
      )
      abline(v = c(-EpsilonR, EpsilonR), lty = 2)
    }

    par(ask = FALSE)
  }

  # Show a summary of the results
  .ShowTestResults(
    ResultDiffMean, EpsilonM, OptThresholdM,
    GroupLabel1, GroupLabel2, "mean expression"
  )
  .ShowTestResults(ResultDiffDisp, EpsilonD, OptThresholdD,
    GroupLabel1, GroupLabel2,
    Task = "over dispersion", Others = NotDE
  )
  if (!is.null(Chain1@parameters$epsilon)) {
    .ShowTestResults(ResultDiffResDisp, EpsilonR, OptThresholdE,
      GroupLabel1, GroupLabel2,
      Task = "residual over dispersion", Others = NotExcluded
    )
  }

  # Output list
  out <- list(
    TableMean = TableMean,
    TableDisp = TableDisp,
    DiffMeanSummary = list(
      ProbThreshold = OptThresholdM[1],
      EFDR = OptThresholdM[2],
      EFNR = OptThresholdM[3]
    ),
    DiffDispSummary = list(
      ProbThreshold = OptThresholdD[1],
      EFDR = OptThresholdD[2],
      EFNR = OptThresholdD[3]
    ),
    Chain1_offset = Chain1_offset,
    Chain2_offset = Chain2_offset,
    OffsetChain = OffsetChain,
    Offset = OffsetEst
  )

  if (!is.null(Chain1@parameters$epsilon)) {
    out$TableResDisp <- TableResDisp
    out$DiffResDispSummary <- list(
      ProbThreshold = OptThresholdE[1],
      EFDR = OptThresholdE[2],
      EFNR = OptThresholdE[3]
    )
  }

  return(out)
}
