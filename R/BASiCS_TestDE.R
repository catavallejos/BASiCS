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
#' @param PlotOffsetType See argument \code{Type} in 
#' \code{\link{BASiCS_PlotOffset}} for more information.
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
#' @param MinESS The minimum effective sample size for a gene to be included 
#' in the tests for
#' differential expression. This helps to remove genes with poor mixing from
#' differential expression tests.
#' Default is 100. If set to NA, genes are
#' not checked for effective sample size before differential expression tests
#' are performed.
#' @param ... Optional parameters.
#'
#' @return \code{BASiCS_TestDE} returns an object of class
#' \code{\linkS4class{BASiCS_ResultsDE}}
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
#' head(format(Test, Which = "Mean"))
#'
#' # Results for the differential over-dispersion test
#' # This only includes genes marked as 'NoDiff' in Test$TableMean
#' head(format(Test, Which = "Disp"))
#'
#' # For testing differences in residual over-dispersion, two chains obtained
#' # via 'BASiCS_MCMC(Data, N, Thin, Burn, Regression=TRUE)' need to be provided
#' data(ChainSCReg)
#' data(ChainRNAReg)
#'
#' Test <- BASiCS_TestDE(
#'   Chain1 = ChainSCReg, Chain2 = ChainRNAReg,
#'   GroupLabel1 = 'SC', GroupLabel2 = 'P&S',
#'   EpsilonM = log2(1.5), EpsilonD = log2(1.5),
#'   EpsilonR = log2(1.5)/log2(exp(1)),
#'   OffSet = TRUE
#' )
#'
#' ## Plotting the results of these tests
#' BASiCS_PlotDE(Test)
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
                          OrderVariable = c("GeneIndex", "GeneName", "Mu"),
                          GroupLabel1 = "Group1",
                          GroupLabel2 = "Group2",
                          Plot = TRUE,
                          PlotOffset = TRUE,
                          PlotOffsetType = c(
                            "offset estimate", 
                            "before-after",
                            "MA plot"
                          ),
                          Offset = TRUE,
                          EFDR_M = 0.05,
                          EFDR_D = 0.05,
                          EFDR_R = 0.05,
                          GenesSelect = rep(
                            TRUE,
                            ncol(Chain1@parameters[["mu"]])
                          ),
                          min.mean = 1,
                          MinESS = 100,
                          ...) {
  
  OrderVariable <- match.arg(OrderVariable)

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
    Offset = Offset
  )

  GeneName <- colnames(Chain1@parameters$mu)
  GeneIndex <- seq_len(length(GeneName))

  message(
    "-------------------------------------------------------------\n",
    "Log-fold change thresholds are now set in a log2 scale. \n",
    "Original BASiCS release used a natural logarithm scale."
  )

  if (xor(
        is.null(Chain1@parameters[["epsilon"]]),
        is.null(Chain2@parameters[["epsilon"]]))
      ) {

    stop("Both chains should be run with the same setting for Regression.")
  }

  n1 <- ncol(Chain1@parameters$nu)
  n2 <- ncol(Chain2@parameters$nu)
  n <- n1 + n2

  IncludeEpsilon <- !is.null(Chain1@parameters[["epsilon"]])
  # With offset correction
  if (Offset) {
    A <- BASiCS_CorrectOffset(Chain1, Chain2, min.mean = min.mean)
    OffsetEst <- A$Offset
    Chain1_offset <- A$Chain
    Chain2_offset <- Chain2  # Chain2 requires no change

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

  } else {
    message(
      "-------------------------------------------------------------\n",
      "It is recomended to perform a global offset correction \n",
      "to remove global changes between the two groups of cells \n",
      "Default offset value set equal to 1.\n",
      "To perform offset correction, please set 'Offset = TRUE'. \n",
      "-------------------------------------------------------------\n"
    )
    Chain1_offset <- Chain1
    Chain2_offset <- Chain2

    # Default values when no offset correction is applied
    OffsetEst <- 1
  }

  Mu1 <- matrixStats::colMedians(Chain1_offset@parameters$mu)
  Mu2 <- matrixStats::colMedians(Chain2_offset@parameters$mu)
  Delta1 <- matrixStats::colMedians(Chain1_offset@parameters$delta)
  Delta2 <- matrixStats::colMedians(Chain2_offset@parameters$delta)
  ChainTau <- log2(Chain1_offset@parameters$mu / Chain2_offset@parameters$mu)

  MuBase <- (Mu1 * n1 + Mu2 * n2) / n
  
  orderVar <- switch(
    OrderVariable,
    "GeneIndex" = order(GeneIndex, decreasing = FALSE),
    "GeneName" = order(GeneName, decreasing = TRUE),
    "Mu" = order(as.numeric(MuBase), decreasing = TRUE)
  )
  q <- length(GenesSelect)

  GoodESS <- .CheckESS(Chain1, Chain2, MinESS, "mu", q)
  ResM <- .RunTest(
    Chain = ChainTau,
    Epsilon = EpsilonM,
    ProbThreshold = ProbThresholdM,
    EFDR = EFDR_M,
    Task = "Differential mean",
    Suffix = "M",
    GroupLabel1 = GroupLabel1,
    GroupLabel2 = GroupLabel2,
    GenesSelect = GenesSelect,
    Param1 = Mu1,
    Param2 = Mu2,
    n1 = n1,
    n2 = n2,
    GeneName = GeneName,
    GoodESS = GoodESS,
    Measure = "Mean"
  )

  # Genes with no change in mean expression
  DE <- !(ResM@Table$ResultDiffMean %in% c("NoDiff", "ExcludedByUser"))

  GoodESS <- .CheckESS(Chain1, Chain2, MinESS, "delta", q)

  ChainOmega <- log2(Chain1@parameters$delta / Chain2@parameters$delta)
  ResD <- .RunTest(
    Chain = ChainOmega,
    Epsilon = EpsilonD,
    ProbThreshold = ProbThresholdD,
    EFDR = EFDR_D,
    Task = "Differential dispersion",
    Suffix = "D",
    GroupLabel1 = GroupLabel1,
    GroupLabel2 = GroupLabel2,
    GenesSelect = GenesSelect,
    Param1 = Delta1,
    Param2 = Delta2,
    n1 = n1,
    n2 = n2,
    GeneName = GeneName,
    Measure = "Disp",
    GoodESS = GoodESS,
    Excluded = DE
  )

  # Changes in residual over-dispersion - if regression approach was used
  if (IncludeEpsilon) {

    Excluded <- is.na(Chain1@parameters$epsilon[1, ]) |
      is.na(Chain2@parameters$epsilon[1, ])

    ChainPsi <- Chain1@parameters$epsilon - Chain2@parameters$epsilon
    Epsilon1 <- matrixStats::colMedians(Chain1@parameters$epsilon)
    Epsilon2 <- matrixStats::colMedians(Chain2@parameters$epsilon)

    GoodESS <- .CheckESS(Chain1, Chain2, MinESS, "epsilon", q)
    ResR <- .RunTest(
      Chain = ChainPsi,
      Epsilon = EpsilonR,
      ProbThreshold = ProbThresholdR,
      EFDR = EFDR_R,
      Task = "Differential residual dispersion",
      Suffix = "R",
      GroupLabel1 = GroupLabel1,
      GroupLabel2 = GroupLabel2,
      GenesSelect = GenesSelect,
      Param1 = Epsilon1,
      Param2 = Epsilon2,
      n1 = n1,
      n2 = n2,
      GeneName = GeneName,
      Measure = "ResDisp",
      GoodESS = GoodESS,
      Excluded = Excluded
    )
    ResR@Table[, "MeanOverall"] <- ResM@Table[, "MeanOverall"]
    ResR <- ResR[orderVar, ]
  }
  ResM <- ResM[orderVar, ]
  ResD <- ResD[orderVar, ]

  Results <- list(
    Mean = ResM,
    Disp = ResD
  )

  if (IncludeEpsilon) {
    Results <- c(
      Results, 
      ResDisp = ResR
    )
  }

  Out <- new("BASiCS_ResultsDE",
    Results = Results,
    Chain1 = Chain1_offset,
    Chain2 = Chain2_offset,
    GroupLabel1 = GroupLabel1,
    GroupLabel2 = GroupLabel2,
    Offset = OffsetEst,
    RowData = DataFrame(GeneName = GeneName),
    Extras = list()
  )
  if (Plot) {
    Out@Extras[["Plots"]] <- BASiCS_PlotDE(Out)
  }
  Out
}


.DiffRes <- function(ResultDE) {
  ResultDE[.WhichDiffRes(ResultDE), ]
}

.WhichDiffRes <- function(ResultDE) {
  !ResultDE@Table$Result %in% c(
    "NoDiff",
    "ExcludedLowESS",
    "ExcludedFromTesting",
    "ExcludedByUser"
  )
}
