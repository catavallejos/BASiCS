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

.ThresholdSearch <- function(Probs, ProbThreshold, EFDR, Task, Suffix)
{
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
    
    EFDRgrid <- vapply(ProbThresholds, FUN = .EFDR, 
                       FUN.VALUE = 1, Prob = Probs)
    EFNRgrid <- vapply(ProbThresholds, FUN = .EFNR, 
                       FUN.VALUE = 1, Prob = Probs)
    
    AbsDiff <- abs(EFDRgrid - EFDR)
    
    if (sum(!is.na(AbsDiff)) > 0) {
      # 1.1. If the calibration doesn't completely fail
      
      # Search EFDR closest to the desired value
      EFDRopt <- EFDRgrid[AbsDiff  == min(AbsDiff , na.rm = TRUE) & !is.na(AbsDiff)]
      # If multiple threholds lead to same EFDR, choose the one with lowest EFNR
      EFNRopt <- EFNRgrid[EFDRgrid == mean(EFDRopt) & !is.na(EFDRgrid)]
      if (sum(!is.na(EFNRopt)) > 0) {
        optimal <- which(EFDRgrid == mean(EFDRopt) & EFNRgrid == mean(EFNRopt))
      } else { optimal <- which(EFDRgrid == mean(EFDRopt)) }
      # Quick fix for EFDR/EFNR ties; possibly not an issue in real datasets
      optimal <- median(round(median(optimal)))
      
      if(ProbThresholds[optimal] >= ProbThreshold) {
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
        message("For ", Task, " task:\n",
                "the posterior probability threshold chosen via EFDR calibration is too low.",
                "Probability threshold automatically set equal to 'ProbThreshold", Suffix, "'.")
      }
    } else {
      # 1.2 If calibration completely fails - default prob to 0.9 (conservative)
      message("EFDR calibration failed for ", Task, " task. \n",
              "Probability threshold automatically set equal to 0.90 \n")
      OptThreshold <- c(0.9, NA, NA)
    }
  } else {
    # 2. If EFDR is not provided - fix to given probs
    
    EFDRgrid <- .EFDR(ProbThreshold, Probs)
    EFNRgrid <- .EFNR(ProbThreshold, Probs)
    
    OptThreshold <- c(ProbThreshold, EFDRgrid[1], EFNRgrid[1])
  }
  
  list(OptThreshold = OptThreshold,
       EFDRgrid = EFDRgrid, 
       EFNRgrid = EFNRgrid)
}
