# Used within BASiCS_TestDE function
HiddenThresholdSearchTestDE <- function(ChainLFC, 
                                        Epsilon,
                                        EviThreshold,
                                        GenesSelect,
                                        EFDR,
                                        Task)
{
  # Calculating posterior probabilities
  if(Epsilon > 0) {Prob = HiddenProbDE(chain = abs(ChainLFC), tol = Epsilon)}
  else 
  {
    Prob_aux = HiddenProbDE(chain = ChainLFC, tol = 0)
    Prob = 2*pmax(Prob_aux, 1-Prob_aux) - 1
  }
  
  # Posterior probability threshold search
  if(is.null(EviThreshold))
  {
    EviThresholds <- seq(0.5,0.9995,by=0.00025)
    
    if(is.null(GenesSelect))
    {
      EFDRgrid <- sapply(EviThresholds, HiddenEFDRDV, Prob = Prob)
      EFNRgrid <- sapply(EviThresholds, HiddenEFNRDV, Prob = Prob)      
    }
    else
    {
      EFDRgrid <- sapply(EviThresholds, HiddenEFDRDV, Prob = Prob[GenesSelect])
      EFNRgrid <- sapply(EviThresholds, HiddenEFNRDV, Prob = Prob[GenesSelect])        
    }
    
    optimal = round(median(which(abs(EFDRgrid - EFDR) == min(abs(EFDRgrid - EFDR)))))
    OptThreshold <- c(EviThresholds[optimal], EFDRgrid[optimal], EFNRgrid[optimal])
    
    EviThreshold = OptThreshold[1]
    if(is.na(EviThreshold)) 
    { 
      message(paste(Task, ": EFDR calibration failed. Probability threshold automatically set equal to 0.90 \n"))
      EviThreshold = 0.90
    }
  }
  else
  {
    if(is.null(GenesSelect))
    {
      EFDRgrid = HiddenEFDRDV(EviThreshold, Prob)
      EFNRgrid = HiddenEFNRDV(EviThreshold, Prob)
    }
    else
    {
      EFDRgrid = HiddenEFDRDV(EviThreshold, Prob[GenesSelect])
      EFNRgrid = HiddenEFNRDV(EviThreshold, Prob[GenesSelect])      
    }
    OptThreshold <- c(EviThreshold, EFDRgrid[1], EFNRgrid)
  }  
  list("Prob" = Prob,
       "OptThreshold" = OptThreshold,
       "EFDRgrid" = EFDRgrid,
       "EFNRgrid" = EFNRgrid)
}