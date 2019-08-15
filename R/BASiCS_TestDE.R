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
#' Default: \code{GenesSelect = NULL}. When used, this argument must be a vector
#' of \code{TRUE} (include gene) / \code{FALSE} (exclude gene) indicator,
#' with the same length as the number of intrinsic genes and following the same
#' order as how genes are displayed in the table of counts.
#' This argument is necessary in order to have a meaningful EFDR calibration
#' when the user decides to exclude some genes from the comparison.
#' @param IncludeEpsilon Flag indicating where to include the residual 
#' overdispersion parameter, epsilon, in the 
#' figures and analysis. Default is \code{TRUE} if epsilon is present in
#' the chains.
#' @param ... Graphical parameters (see \code{\link[graphics]{par}}).
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
#' Test <- BASiCS_TestDE(Chain1 = ChainSC, Chain2 = ChainRNA,
#'                       GroupLabel1 = 'SC', GroupLabel2 = 'P&S',
#'                       EpsilonM = log2(1.5), EpsilonD = log2(1.5),
#'                       OffSet = TRUE)
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
#' Test <- BASiCS_TestDE(Chain1 = ChainSCReg, Chain2 = ChainRNAReg,
#'                       GroupLabel1 = 'SC', GroupLabel2 = 'P&S',
#'                       EpsilonM = log2(1.5), EpsilonD = log2(1.5),
#'                       EpsilonR = log2(1.5)/log2(exp(1)),
#'                       OffSet = TRUE)
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
                          GenesSelect = NULL, 
                          IncludeEpsilon = !is.null(Chain1@parameters[["epsilon"]]),
                          ...)
{
  OrderVariable <- match.arg(OrderVariable)

  HiddenHeaderTest_DE(Chain1,
                      Chain2,
                      EpsilonM,
                      EpsilonD,
                      EpsilonR,
                      EFDR_M,
                      EFDR_D,
                      EFDR_R,
                      ProbThresholdM,
                      ProbThresholdD,
                      ProbThresholdR,
                      OrderVariable,
                      GroupLabel1,
                      GroupLabel2,
                      GenesSelect,
                      Plot,
                      PlotOffset,
                      Offset)

  GeneName <- colnames(Chain1@parameters$mu)
  GeneIndex <- seq_len(length(GeneName))

  message("-------------------------------------------------------------\n",
          "Log-fold change thresholds are now set in a log2 scale. \n",
          "Original BASiCS release used a natural logarithm scale.")

  n1 <- ncol(Chain1@parameters$nu)
  n2 <- ncol(Chain2@parameters$nu)
  n <- n1 + n2

  # With offset correction
  if (Offset) {

    OffsetCorrected <- BASiCS_CorrectOffset(Chain1 = Chain1, 
                                            Chain2 = Chain2, 
                                            GroupLabel1 = GroupLabel1, 
                                            GroupLabel2 = GroupLabel2, 
                                            Plot = PlotOffset)
    Mu1 <- OffsetCorrected@Mu1
    Mu2 <- OffsetCorrected@Mu2
    Delta1 <- OffsetCorrected@Delta1
    Delta2 <- OffsetCorrected@Delta2
    MuBase <- OffsetCorrected@MuBase
    ChainTau <- OffsetCorrected@ChainTau
    MedianTau <- OffsetCorrected@MedianTau

    # Default values when no offset correction is applied
    OffsetEst <- OffsetCorrected@OffsetEst
    OffsetChain <- OffsetCorrected@OffsetChain
    Chain1_offset <- OffsetCorrected@Chain1
    Chain2_offset <- OffsetCorrected@Chain2

  } else {
    OffsetChain <- matrixStats::rowSums2(Chain1@parameters$mu) /
                   matrixStats::rowSums2(Chain2@parameters$mu)
    # Offset point estimate
    OffsetEst <- median(OffsetChain)
    if (!isTRUE(all.equal(OffsetEst, 0))) {
      stop(paste("Global offset detected between Chain1 and Chain2!\n ",
                 "Please remove with BASiCS_correctOffset or set Offset=TRUE."))
    }

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



  TestDifferential <- function(Chain, 
                               Epsilon,
                               Base,
                               Param1,
                               Param2,
                               ProbThreshold, 
                               GenesSelect, 
                               EFDR,
                               GroupLabel1,
                               GroupLabel2, 
                               Aux,
                               Measure
                               ) {

    Median <- matrixStats::colMedians(Chain)
    Base <- (Param1 * n1 + Param2 * n2) / n

    Prob <- Aux$Prob
    OptThreshold <- Aux$OptThreshold

    # Test results
    Plus1 <- which(Prob > OptThreshold[[1]] & Median > 0)
    Plus2 <- which(Prob > OptThreshold[[1]] & Median < 0)
    ResultDiff <- rep("NoDiff", length(Median))
    ResultDiff[Plus1] <- paste0(GroupLabel1, "+")
    ResultDiff[Plus2] <- paste0(GroupLabel2, "+")
    if (!is.null(GenesSelect)) {
        ResultDiff[!GenesSelect] <- "ExcludedByUser"
    }
    # Output table
    Table <- cbind.data.frame(
      GeneName = GeneName,
      MEASUREOverall = as.numeric(Base),
      MEASURE1 = Param1,
      MEASURE2 = Param2,
      MEASUREFC = as.numeric(2 ^ Median),
      MEASUREDISTANCE = as.numeric(Median),
      ProbDiffMEASURE = as.numeric(Prob),
      ResultDiffMEASURE = ResultDiff,
      stringsAsFactors = FALSE)

    if (Measure == "epsilon") {
      Table$MeasureFC <- NULL
    }
    colnames(Table) <- gsub("MEASURE", Measure, colnames(Table))
    colnames(Table) <- gsub("DISTANCE", DistanceVar(Measure), colnames(Table))

    # Rounding to 3 decimal points
    IndNumeric <- vapply(Table, is.numeric, logical(1))
    Table[, IndNumeric] <- round(Table[, IndNumeric], 3)
    
    # PlotSearch(Measure, Aux, EFDR)

    Table
  }


  AuxMean <- HiddenThresholdSearchTestDE(ChainTau,
                                         EpsilonM,
                                         ProbThresholdM,
                                         GenesSelect,
                                         EFDR_M,
                                         Task = "Differential expression", 
                                         Suffix = "M")

  TableMean <- TestDifferential(
    Chain = ChainTau, 
    Epsilon = EpsilonM, 
    Param1 = Mu1,
    Param2 = Mu2,
    ProbThreshold = ProbThresholdM, 
    GenesSelect = GenesSelect, 
    EFDR = EFDR_M, 
    GroupLabel1 = GroupLabel1,
    GroupLabel2 = GroupLabel2,
    Aux = AuxMean,
    Measure = "Mean"
  )

  # Genes with no change in mean expression
  NotDE <- TableMean$ResultDiffMean == "NoDiff"

  # Genes to calibrate EFDR
  if (!is.null(GenesSelect)) {
    select <- NotDE & GenesSelect
  } else {
    select <- NotDE
  }

  ChainOmega <- log2(Chain1@parameters$delta / Chain2@parameters$delta)


  AuxDisp <- HiddenThresholdSearchTestDE(
    ChainOmega,
    EpsilonD,
    ProbThresholdD,
    select,
    EFDR_D,
    Task = "Differential dispersion", 
    Suffix = "D")


  TableDisp <- TestDifferential(
    Chain = ChainOmega, 
    Epsilon = EpsilonD, 
    Param1 = Delta1,
    Param2 = Delta2,
    ProbThreshold = ProbThresholdD, 
    GenesSelect = select, 
    EFDR = EFDR_D, 
    GroupLabel1 = GroupLabel1,
    GroupLabel2 = GroupLabel2,
    Aux = AuxDisp,
    Measure = "Disp"
  )

  TableDisp$MeanOverall <- TableMean$MeanOverall
  TableDisp$ResultDiffDisp[!NotDE] <- "ExcludedFromTesting"
  if (!is.null(GenesSelect)) {
    TableDisp$ResultDiffDisp[!GenesSelect] <- "ExcludedByUser"
  }

  # Changes in over dispersion
  orderVar <- switch(OrderVariable,
    "GeneIndex" = order(GeneIndex, decreasing = FALSE),
    "GeneName" = order(GeneName, decreasing = TRUE),
    "Mu" = order(as.numeric(MuBase), decreasing = TRUE)
  )

  # Changes in residual over-dispersion - if regression approach was used
  if (IncludeEpsilon) {

    NotExcluded <- !(is.na(Chain1@parameters$epsilon[1, ]) |
                     is.na(Chain2@parameters$epsilon[1, ]))
    # Genes to calibrate EFDR
    if (!is.null(GenesSelect)) {
      select <- NotExcluded & GenesSelect
    } else {
      select <- NotExcluded
    }

    ChainPsi <- Chain1@parameters$epsilon - Chain2@parameters$epsilon
    Epsilon1 <- matrixStats::colMedians(Chain1@parameters$epsilon)
    Epsilon2 <- matrixStats::colMedians(Chain2@parameters$epsilon)

    AuxResDisp <- HiddenThresholdSearchTestDE(
      ChainPsi,
      EpsilonR,
      ProbThresholdR,
      select,
      EFDR_R,
      Task = "residual differential dispersion", 
      Suffix = "R")


    TableResDisp <- TestDifferential(
      Chain = ChainPsi, 
      Epsilon = EpsilonR, 
      Param1 = Epsilon1,
      Param2 = Epsilon2,
      ProbThreshold = ProbThresholdR, 
      GenesSelect = select,
      EFDR = EFDR_R, 
      GroupLabel1 = GroupLabel1,
      GroupLabel2 = GroupLabel2,
      Aux = AuxResDisp,
      Measure = "ResDisp")

    TableResDisp$MeanOverall <- TableMean$MeanOverall
    TableResDisp$ResultDiffResDisp[!NotExcluded] <- "ExcludedFromTesting"
    if (!is.null(GenesSelect)) {
      TableResDisp$ResultDiffResDisp[!GenesSelect] <- "ExcludedByUser"
    }

    TableResDisp <- TableResDisp[orderVar, ]
  }

  TableMean <- TableMean[orderVar, ]
  TableDisp <- TableDisp[orderVar, ]

  if (!is.null(GenesSelect)) {
    message("-------------------------------------------------------------\n",
            "The user excluded ", sum(!GenesSelect), " genes. \n",
            "These genes are marked as 'ExcludedByUser' \n",
            "and excluded from EFDR calibration. \n",
            "-------------------------------------------------------------\n")
  }

  Results <- list(
    Mean = new("BASiCS_ResultDE", 
      Table = TableMean,
      Name = "Mean",
      GroupLabel1 = GroupLabel1,
      GroupLabel2 = GroupLabel2,
      ProbThreshold = AuxMean$OptThreshold[[1]],
      EFDR = AuxMean$OptThreshold[[2]],
      EFNR = AuxMean$OptThreshold[[3]],
      EFDRgrid = AuxMean$EFDRgrid,
      EFNRgrid = AuxMean$EFNRgrid,
      Epsilon = EpsilonM
    ),
    Disp = new("BASiCS_ResultDE", 
      Table = TableDisp,
      Name = "Disp",
      GroupLabel1 = GroupLabel1,
      GroupLabel2 = GroupLabel2,
      ProbThreshold = AuxDisp$OptThreshold[[1]],
      EFDR = AuxDisp$OptThreshold[[2]],
      EFNR = AuxDisp$OptThreshold[[3]],
      EFDRgrid = AuxDisp$EFDRgrid,
      EFNRgrid = AuxDisp$EFNRgrid,
      Epsilon = EpsilonD
    )
  )

  if (IncludeEpsilon) {
    Results <- c(
      Results, 
      ResDisp = new("BASiCS_ResultDE", 
        Table = TableResDisp,
        Name = "ResDisp",
        GroupLabel1 = GroupLabel1,
        GroupLabel2 = GroupLabel2,
        ProbThreshold = AuxResDisp$OptThreshold[[1]],
        EFDR = AuxResDisp$OptThreshold[[2]],
        EFNR = AuxResDisp$OptThreshold[[3]],
        EFDRgrid = AuxResDisp$EFDRgrid,
        EFNRgrid = AuxResDisp$EFNRgrid,
        Epsilon = EpsilonR
      )
    )
  }



  Out <- new("BASiCS_ResultsDE",
    Results = Results,
    Chain1 = Chain1_offset,
    Chain2 = Chain2_offset,
    GroupLabel1 = GroupLabel1,
    GroupLabel2 = GroupLabel2,
    OffsetChain = OffsetChain,
    Offset = OffsetEst,
    Extras = list()
  )
  if (Plot) {
    BASiCS_PlotDE(Out)
  }
  Out
}

