.TailProb <- function(Chain, Threshold) {
  if (Threshold > 0) {
    Prob <- matrixStats::colMeans2(abs(Chain) > Threshold)
  }
  else {
    Prob_aux <- matrixStats::colMeans2(Chain > 0)
    Prob <- 2 * pmax(Prob_aux, 1 - Prob_aux) - 1
  }
  return(Prob)
}

.EFDR <- function(EviThreshold, Prob) {
  return(sum((1 - Prob) * I(Prob > EviThreshold)) / sum(I(Prob > EviThreshold)))
}

.EFNR <- function(EviThreshold, Prob) {
  return(sum(Prob * I(EviThreshold >= Prob)) / sum(I(EviThreshold >= Prob)))
}

.ThresholdSearch <- function(Probs, ProbThreshold, EFDR, Task, Suffix) {
  # Summary of cases
  # 1. If EFDR is provided - run calibration
  #   1.1. If the calibration doesn't completely fail - search prob
  #     1.1.1. If optimal prob is not too low - set prob to optimal
  #     1.1.2. If optimal prob is too low - fix to input probs
  #   1.2 If calibration completely fails - default prob to 0.9 (conservative)
  # 2. If EFDR is not provided - fix to input probs
  
  # 2. If EFDR is not provided - fix to given probs
  if (is.null(EFDR)) {
    EFDRgrid <- .EFDR(ProbThreshold, Probs)
    EFNRgrid <- .EFNR(ProbThreshold, Probs)
    
    OptThreshold <- c(ProbThreshold, EFDRgrid[1], EFNRgrid[1])
    return(
      list(
        OptThreshold = OptThreshold,
        EFDRgrid = EFDRgrid, 
        EFNRgrid = EFNRgrid
      )
    )
  }



# 1. If EFDR is provided - run calibration    
  ProbThresholds <- seq(0.5, 0.9995, by = 0.00025)
  
  EFDRgrid <- vapply(
    ProbThresholds,
    FUN = .EFDR, 
    FUN.VALUE = 1,
    Prob = Probs
  )
  EFNRgrid <- vapply(
    ProbThresholds,
    FUN = .EFNR,
    FUN.VALUE = 1,
    Prob = Probs
  )
  
  AbsDiff <- abs(EFDRgrid - EFDR)
  

  # 1.2 If calibration completely fails - default prob to 0.9 (conservative)
  if (sum(!is.na(AbsDiff)) == 0) {
    message(
      "EFDR calibration failed for ", Task, " task. \n",
      "Probability threshold automatically set equal to 0.90 \n"
    )
    OptThreshold <- c(0.9, NA, NA)
    return(
      list(
        OptThreshold = OptThreshold,
        EFDRgrid = EFDRgrid, 
        EFNRgrid = EFNRgrid
      )
    )
  }


  # 1.1. If the calibration doesn't completely fail
  
  # Search EFDR closest to the desired value
  EFDRopt <- EFDRgrid[AbsDiff  == min(AbsDiff , na.rm = TRUE) & 
    !is.na(AbsDiff)]
  # If multiple threholds lead to same EFDR, choose the one with lowest EFNR
  EFNRopt <- EFNRgrid[EFDRgrid == mean(EFDRopt) & !is.na(EFDRgrid)]
  if (sum(!is.na(EFNRopt)) > 0) {
    optimal <- which(EFDRgrid == mean(EFDRopt) & EFNRgrid == mean(EFNRopt))
  } else {
    optimal <- which(EFDRgrid == mean(EFDRopt))
  }
  # Quick fix for EFDR/EFNR ties; possibly not an issue in real datasets
  optimal <- median(round(median(optimal)))
  
  if (ProbThresholds[optimal] >= ProbThreshold) {
    # 1.1.1. If optimal prob is not too low - set prob to optimal
    
    OptThreshold <- c(ProbThresholds[optimal],
                      EFDRgrid[optimal], 
                      EFNRgrid[optimal])
    if (abs(OptThreshold[2] - EFDR) > 0.025) {
      message("For ", Task, " task:\n",
              "It is not possible to find a probability threshold (>0.5) \n",
              "that achieves the desired EFDR level (+-0.025). \n",
              "The output below reflects the closest possible value. \n")
    }
  } else {
    # 1.1.2. If optimal prob is too low - fix to input probs
    
    EFDRgrid <- .EFDR(ProbThreshold, Probs)
    EFNRgrid <- .EFNR(ProbThreshold, Probs)
    
    OptThreshold <- c(ProbThreshold, EFDRgrid[1], EFNRgrid[1])  
    message(
      "For ", Task, " task:\n",
      "the posterior probability threshold chosen via EFDR calibration",
      "is too low. Probability threshold automatically set equal to",
      "'ProbThreshold", Suffix, "'.")
  }

  list(
    OptThreshold = OptThreshold,
    EFDRgrid = EFDRgrid, 
    EFNRgrid = EFNRgrid
  )
}

.TestDifferential <- function(
    Chain,
    Param1,
    Param2,
    n1,
    n2,
    GenesSelect,
    GroupLabel1,
    GroupLabel2, 
    Prob,
    OptThreshold,
    GeneName,
    Measure,
    GoodESS,
    Excluded = NULL
  ) {

  Median <- matrixStats::colMedians(Chain)
  Base <- (Param1 * n1 + Param2 * n2) / (n1 + n2)

  Result <- .TestResults(
    Prob = Prob,
    Threshold = OptThreshold[[1]],
    Estimate = Median,
    Label1 = GroupLabel1,
    Label2 = GroupLabel2,
    GenesSelect = GenesSelect,
    GoodESS = GoodESS, 
    Excluded = Excluded
  )
  # Output table
  Table <- cbind.data.frame(
    GeneName = GeneName,
    MEASUREOverall = as.numeric(Base),
    MEASURE1 = Param1,
    MEASURE2 = Param2,
    MEASUREFC = as.numeric(2 ^ Median),
    MEASUREDISTANCE = as.numeric(Median),
    ProbDiffMEASURE = as.numeric(Prob),
    ResultDiffMEASURE = Result,
    stringsAsFactors = FALSE
  )

  if (Measure == "epsilon") {
    Table$MeasureFC <- NULL
  }
  colnames(Table) <- gsub("MEASURE", Measure, colnames(Table))
  colnames(Table) <- gsub("DISTANCE", .DistanceVar(Measure), colnames(Table))

  # Rounding to 3 decimal points
  IndNumeric <- vapply(Table, is.numeric, logical(1))
  Table[, IndNumeric] <- round(Table[, IndNumeric], digits = 3)
  
  Table
}


.TestResults <- function(Prob,
                         Threshold,
                         Estimate,
                         Label1,
                         Label2,
                         GenesSelect,
                         GoodESS,
                         Excluded = NULL) {

  # Which genes are + in each group
  Plus1 <- which(Prob > Threshold & Estimate > 0)
  Plus2 <- which(Prob > Threshold & Estimate < 0)
  Result <- rep("NoDiff", length(Estimate))
  Result[Plus1] <- paste0(Label1, "+")
  Result[Plus2] <- paste0(Label2, "+")
  if (!is.null(Excluded)) {
    Result[Excluded] <- "ExcludedFromTesting"
  }
  Result[!GenesSelect] <- "ExcludedByUser"
  Result[!GoodESS] <- "ExcludedLowESS"
  return(Result)
}

.RunTest <- function(Chain,
                     Epsilon,
                     ProbThreshold,
                     EFDR,
                     Task,
                     Suffix,
                     GroupLabel1,
                     GroupLabel2,
                     GenesSelect,
                     Param1,
                     Param2,
                     n1,
                     n2,
                     GeneName,
                     Measure,
                     GoodESS,
                     Excluded = NULL) {

  Median <- matrixStats::colMedians(Chain)

  Prob <- .TailProb(Chain = Chain, Threshold = Epsilon)
  if (!is.null(Excluded)) {
    select <- GenesSelect & GoodESS & !Excluded
  } else {
    select <- GenesSelect & GoodESS
  }

  Aux <- .ThresholdSearch(
    Probs = Prob[select],
    ProbThreshold = ProbThreshold,
    EFDR = EFDR,
    Task = Task,
    Suffix = Suffix
  )
  OptThreshold <- Aux$OptThreshold

  Table <- .TestDifferential(
    Chain = Chain,
    Param1 = Param1,
    Param2 = Param2,
    n1 = n1,
    n2 = n2,
    GenesSelect = GenesSelect, 
    GroupLabel1 = GroupLabel1,
    GroupLabel2 = GroupLabel2,
    Prob = Prob,
    OptThreshold = OptThreshold,
    GeneName = GeneName,
    Measure = Measure,
    GoodESS = GoodESS,
    Excluded = Excluded
  )

  new(
    "BASiCS_ResultDE", 
    Table = Table,
    Name = Measure,
    GroupLabel1 = GroupLabel1,
    GroupLabel2 = GroupLabel2,
    ProbThreshold = Aux$OptThreshold[[1]],
    EFDR = Aux$OptThreshold[[2]],
    EFNR = Aux$OptThreshold[[3]],
    EFDRgrid = Aux$EFDRgrid,
    EFNRgrid = Aux$EFNRgrid,
    Epsilon = Epsilon
  )
}

.CheckESS <- function(Chain1, Chain2, MinESS, parameter, q) {
  if (!is.na(MinESS)) {
    ess(mcmc(Chain1@parameters[[parameter]])) > MinESS &
    ess(mcmc(Chain2@parameters[[parameter]])) > MinESS
  } else {
    rep(TRUE, q)
  }
}
