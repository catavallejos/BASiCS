HiddenThresholdSearchTestDE <- function(ChainLFC, Epsilon, ProbThreshold, 
                                        GenesSelect, EFDR, Task) 
{
  # Calculating posterior probabilities
  if (Epsilon > 0) { Prob <- matrixStats::colMeans2(abs(ChainLFC) > Epsilon) } 
  else 
  {
    Prob_aux <- matrixStats::colMeans2(ChainLFC > 0)
    Prob <- 2 * pmax(Prob_aux, 1 - Prob_aux) - 1
  }
  
  # Posterior probability threshold search
  if (is.null(ProbThreshold)) 
  {
    # When a posterior probability threshold has not been set a priori
    ProbThresholds <- seq(0.5, 0.9995, by = 0.00025)
    
    if (is.null(GenesSelect)) 
    {
      EFDRgrid <- sapply(ProbThresholds, HiddenEFDRDV, Prob = Prob)
      EFNRgrid <- sapply(ProbThresholds, HiddenEFNRDV, Prob = Prob)
    } 
    else 
    {
      EFDRgrid <- sapply(ProbThresholds, HiddenEFDRDV, Prob = Prob[GenesSelect])
      EFNRgrid <- sapply(ProbThresholds, HiddenEFNRDV, Prob = Prob[GenesSelect])
    }
    
    above <- abs(EFDRgrid - EFDR)
    
    if(sum(!is.na(above)) > 0)
    {
      # Search EFDR closest to the desired value
      EFDRopt <- EFDRgrid[above == min(above, na.rm = TRUE) & !is.na(above)] 
      # If multiple threholds lead to same EFDR, choose the one with lowest EFNR
      EFNRopt <- EFNRgrid[EFDRgrid == mean(EFDRopt) & !is.na(EFDRgrid)]
      if(sum(!is.na(EFNRopt)) > 0)
      {
        optimal <- which(EFDRgrid == mean(EFDRopt) & EFNRgrid == mean(EFNRopt))
      } 
      else
      {
        optimal <- which(EFDRgrid == mean(EFDRopt))
      }
      # Quick fix for EFDR/EFNR ties; possibly not an issue in real datasets
      optimal <- median(round(median(optimal)))
      OptThreshold <- c(ProbThresholds[optimal], 
                        EFDRgrid[optimal], EFNRgrid[optimal])
      
      if( abs(OptThreshold[2] - EFDR) > 0.025 )
      {
        message("For ", Task, " task:\n",
                "It is not possible to find a probability threshold (>0.5) \n",
                "that achieves the desired EFDR level (+-0.025). \n",
                "The output below reflects the closest possible value. \n")         
      }  
    }
    else
    {
      message("EFDR calibration failed for ", Task, "task. \n",
              "Probability threshold automatically set equal to 0.90 \n")
      OptThreshold <- c(0.9, NA, NA)  
    }
  } 
  else 
  {
    # When a posterior probability threshold has been set a priori
    if (is.null(GenesSelect)) 
    {
      EFDRgrid <- HiddenEFDRDV(ProbThreshold, Prob)
      EFNRgrid <- HiddenEFNRDV(ProbThreshold, Prob)
    } 
    else 
    {
      EFDRgrid <- HiddenEFDRDV(ProbThreshold, Prob[GenesSelect])
      EFNRgrid <- HiddenEFNRDV(ProbThreshold, Prob[GenesSelect])
    }
    OptThreshold <- c(ProbThreshold, EFDRgrid[1], EFNRgrid)
  }
  
  list(Prob = Prob, OptThreshold = OptThreshold, 
       EFDRgrid = EFDRgrid, EFNRgrid = EFNRgrid)
}
