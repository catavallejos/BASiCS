# Used in BASiCS_DetectHVG and BASiCS_DetectLVG
HiddenThresholdSearchDetectHVG_LVG <- function(ProbThreshold,
                                               Prob,
                                               EFDR)
{
  if (is.null(ProbThreshold)) {
    # If EviThreshold is not set a priori (search)
    ProbThresholds <- seq(0.5,0.9995,by=0.0005)

    EFDRgrid <- vapply(ProbThresholds, FUN = HiddenEFDR, 
                       FUN.VALUE = 1, Prob = Prob)
    EFNRgrid <- vapply(ProbThresholds, FUN = HiddenEFNR, 
                       FUN.VALUE = 1, Prob = Prob)

    above <- abs(EFDRgrid - EFDR)

    if (sum(!is.na(above)) > 0) {
      # Search EFDR closest to the desired value
      EFDRopt <- EFDRgrid[above == min(above, na.rm = TRUE) & !is.na(above)]
      # If multiple threholds lead to same EFDR, choose the lowest EFDR
      EFNRopt <- EFNRgrid[EFDRgrid == mean(EFDRopt) & !is.na(EFDRgrid)]
      if (sum(!is.na(EFNRopt)) > 0) {
        optimal <- which(EFDRgrid == mean(EFDRopt) & EFNRgrid == mean(EFNRopt))
      } else {
        optimal <- which(EFDRgrid == mean(EFDRopt))
      }
      # Quick fix for EFDR/EFNR ties; possibly not an issue in real datasets
      optimal <- median(round(median(optimal)))
      OptThreshold <- c(ProbThresholds[optimal],
                        EFDRgrid[optimal], EFNRgrid[optimal])

      if (OptThreshold[1] < 0.5) {
        message("For the given variance contribution threshold, \n",
                "the probability threshold with the desired EFDR is below 0.5. \n",
                "Probability threshold set at 0.5 (EFDR/EFNR reported)")
        OptThreshold <- c(ProbThresholds[1], EFDRgrid[1], EFNRgrid[1])
      }

      if (abs(OptThreshold[2] - EFDR) > 0.025) {
        # Message when different to desired EFDR is large
        message("For the given variance contribution threshold, \n",
                "it is not possible to find a probability threshold (>0.5) \n",
                "that achieves the desired EFDR level (tolerance +- 0.025). \n",
                "Output based on the closest possible value. \n")
      }
    }
    else {
      message("For the given variance contribution threshold, \n",
              "it is not possible to estimate EFDR. \n",
              "Probability threshold set at 0.5. \n")
      OptThreshold <- c(0.5, NA, NA)
    }
  }
  else {
    # If EviThreshold is set a priori
    EFDR = HiddenEFDR(ProbThreshold, Prob)
    EFNR = HiddenEFNR(ProbThreshold, Prob)
    OptThreshold <- c(ProbThreshold, EFDR, EFNR)
    ProbThresholds <- NULL
    EFDRgrid <- NULL
    EFNRgrid <- NULL
  }
  list("ProbThresholds" = ProbThresholds,
       "EFDRgrid" = EFDRgrid,
       "EFNRgrid" = EFNRgrid,
       "OptThreshold" = OptThreshold)
}
