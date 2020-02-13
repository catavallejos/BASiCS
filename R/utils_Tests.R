.TailProb <- function(Chain, Threshold) {
  if (Threshold > 0) {
    Prob <- matrixStats::colMeans2(ifelse(Chain > Threshold, 1, 0))
  }
  else {
    Prob_aux <- matrixStats::colMeans2(ifelse(Chain > 0, 1, 0))
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

  # 1. If EFDR is provided - run calibration
  if (!is.null(EFDR)) {
    ProbThresholds <- seq(0.5, 0.9995, by = 0.00025)

    EFDRgrid <- vapply(ProbThresholds,
      FUN = .EFDR,
      FUN.VALUE = 1, Prob = Probs
    )
    EFNRgrid <- vapply(ProbThresholds,
      FUN = .EFNR,
      FUN.VALUE = 1, Prob = Probs
    )

    AbsDiff <- abs(EFDRgrid - EFDR)

    if (sum(!is.na(AbsDiff)) > 0) {
      # 1.1. If the calibration doesn't completely fail

      # Search EFDR closest to the desired value
      EFDRopt <- EFDRgrid[AbsDiff == min(AbsDiff, na.rm = TRUE) & !is.na(AbsDiff)]
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

        OptThreshold <- c(
          ProbThresholds[optimal],
          EFDRgrid[optimal],
          EFNRgrid[optimal]
        )
        if (abs(OptThreshold[2] - EFDR) > 0.025) {
          message(
            "For ", Task, " task:\n",
            "It is not possible to find a probability threshold (>0.5) \n",
            "that achieves the desired EFDR level (+-0.025). \n",
            "The output below reflects the closest possible value. \n"
          )
        }
      } else {
        # 1.1.2. If optimal prob is too low - fix to input probs

        EFDRgrid <- .EFDR(ProbThreshold, Probs)
        EFNRgrid <- .EFNR(ProbThreshold, Probs)

        OptThreshold <- c(ProbThreshold, EFDRgrid[1], EFNRgrid[1])
        message(
          "For ", Task, " task:\n",
          "the posterior probability threshold chosen via EFDR calibration is too low.",
          "Probability threshold automatically set equal to 'ProbThreshold", Suffix, "'."
        )
      }
    } else {
      # 1.2 If calibration completely fails - default prob to 0.9 (conservative)
      message(
        "EFDR calibration failed for ", Task, " task. \n",
        "Probability threshold automatically set equal to 0.90 \n"
      )
      OptThreshold <- c(0.9, NA, NA)
    }
  } else {
    # 2. If EFDR is not provided - fix to given probs

    EFDRgrid <- .EFDR(ProbThreshold, Probs)
    EFNRgrid <- .EFNR(ProbThreshold, Probs)

    OptThreshold <- c(ProbThreshold, EFDRgrid[1], EFNRgrid[1])
  }

  list(
    OptThreshold = OptThreshold,
    EFDRgrid = EFDRgrid,
    EFNRgrid = EFNRgrid
  )
}

.TestResults <- function(Probs,
                         Threshold,
                         Estimate,
                         Label1,
                         Label2,
                         GenesSelect,
                         GoodESS,
                         Excluded = NULL) {

  # Which genes are + in each group
  Plus1 <- which(Probs > Threshold & Estimate > 0)
  Plus2 <- which(Probs > Threshold & Estimate < 0)
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

.ShowTestResults <- function(Result, Epsilon, Opt, Label1, Label2,
                             Task, Others = NULL) {
  n1 <- sum(Result == paste0(Label1, "+"))
  n2 <- sum(Result == paste0(Label2, "+"))

  message("-------------------------------------------------------------\n")

  message(
    n1 + n2, paste0(" genes with a change in ", Task, ":\n"),
    "- Higher in ", Label1, " samples: ", n1, "\n",
    "- Higher in ", Label2, " samples: ", n2, "\n",
    "- Fold change tolerance = ", round(2^(Epsilon) * 100, 2), "% \n",
    "- Probability threshold = ", Opt[1], "\n",
    "- EFDR = ", round(100 * Opt[2], 2), "% \n",
    "- EFNR = ", round(100 * Opt[3], 2), "% \n"
  )

  if (Task == "over dispersion" & !is.null(Others)) {
    message(
      "-------------------------------------------------------------\n",
      "NOTE: differential dispersion assessment only applied to the \n",
      sum(Others), " genes for which the mean did not change. \n",
      "and that were included for testing. \n"
    )
  }
  if (Task == "residual over dispersion" & !is.null(Others)) {
    message(
      "-------------------------------------------------------------\n",
      "NOTE: differential residual dispersion assessment applied to \n",
      sum(Others), " genes expressed in at least 2 cells per condition \n",
      "and that were included for testing. \n"
    )
  }

  message("-------------------------------------------------------------\n\n")
}
