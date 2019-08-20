HiddenThresholdSearchTestDE <- function(ChainLFC, Epsilon, ProbThreshold,
                                        GenesSelect, EFDR, Task, suffix)
{
  # Calculating posterior probabilities
  Prob <- HiddenTailProbTestDE(ChainLFC, Epsilon)

  # Posterior probability threshold search
  if (!is.null(EFDR)) {
    
    ProbThresholds <- seq(0.5, 0.9995, by = 0.00025)

    if (is.null(GenesSelect)) {
      EFDRgrid <- vapply(ProbThresholds, FUN = HiddenEFDR, 
                         FUN.VALUE = 1, Prob = Prob)
      EFNRgrid <- vapply(ProbThresholds, FUN = HiddenEFNR, 
                         FUN.VALUE = 1, Prob = Prob)
    }
    else {
      EFDRgrid <- vapply(ProbThresholds, FUN = HiddenEFDR, 
                         FUN.VALUE = 1, Prob = Prob[GenesSelect])
      EFNRgrid <- vapply(ProbThresholds, FUN = HiddenEFNR, 
                         FUN.VALUE = 1, Prob = Prob[GenesSelect])
    }

    above <- abs(EFDRgrid - EFDR)

    if (sum(!is.na(above)) > 0) {
      
      # Search EFDR closest to the desired value
      EFDRopt <- EFDRgrid[above == min(above, na.rm = TRUE) & !is.na(above)]
      # If multiple threholds lead to same EFDR, choose the one with lowest EFNR
      EFNRopt <- EFNRgrid[EFDRgrid == mean(EFDRopt) & !is.na(EFDRgrid)]
      if (sum(!is.na(EFNRopt)) > 0) {
        optimal <- which(EFDRgrid == mean(EFDRopt) & EFNRgrid == mean(EFNRopt))
      }
      else {
        optimal <- which(EFDRgrid == mean(EFDRopt))
      }
      # Quick fix for EFDR/EFNR ties; possibly not an issue in real datasets
      optimal <- median(round(median(optimal)))
      
      # If calibrated threshold is above the minimum required probability
      if(ProbThresholds[optimal] >= ProbThreshold) {
        
        OptThreshold <- c(ProbThresholds[optimal],
                          EFDRgrid[optimal], EFNRgrid[optimal])
        if (abs(OptThreshold[2] - EFDR) > 0.025) {
          message("For ", Task, " task:\n",
                  "It is not possible to find a probability threshold (>0.5) \n",
                  "that achieves the desired EFDR level (+-0.025). \n",
                  "The output below reflects the closest possible value. \n")
        }
      }
      # If calibrated threshold is below the minimum required probability
      else {
        if (is.null(GenesSelect)) {
          EFDRgrid <- HiddenEFDR(ProbThreshold, Prob)
          EFNRgrid <- HiddenEFNR(ProbThreshold, Prob)
        }
        else {
          EFDRgrid <- HiddenEFDR(ProbThreshold, Prob[GenesSelect])
          EFNRgrid <- HiddenEFNR(ProbThreshold, Prob[GenesSelect])
        }
        OptThreshold <- c(ProbThreshold, EFDRgrid[1], EFNRgrid[1])  
        message("For ", Task, " task:\n",
                "the posterior probability threshold chosen via EFDR calibration is too low.",
                "Probability threshold automatically set equal to 'ProbThreshold", suffix, "'.")
      }
    }
    else {
      message("EFDR calibration failed for ", Task, " task. \n",
              "Probability threshold automatically set equal to 0.90 \n")
      OptThreshold <- c(0.9, NA, NA)
    }
  }
  else {
    # When a posterior probability threshold has been set a priori
    if (is.null(GenesSelect)) {
      EFDRgrid <- HiddenEFDR(ProbThreshold, Prob)
      EFNRgrid <- HiddenEFNR(ProbThreshold, Prob)
    }
    else {
      EFDRgrid <- HiddenEFDR(ProbThreshold, Prob[GenesSelect])
      EFNRgrid <- HiddenEFNR(ProbThreshold, Prob[GenesSelect])
    }
    OptThreshold <- c(ProbThreshold, EFDRgrid[1], EFNRgrid[1])
  }

  list(Prob = Prob, OptThreshold = OptThreshold,
       EFDRgrid = EFDRgrid, EFNRgrid = EFNRgrid)
}
