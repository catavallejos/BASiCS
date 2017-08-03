# Used in BASiCS_DetectHVG and BASiCS_DetectLVG
HiddenThresholdSearchDetectHVG_LVG <- function(EviThreshold,
                                               VarThreshold,
                                               Prob, 
                                               EFDR)
{
  # If EviThreshold is not set a priori (search)
  if(length(EviThreshold) == 0)
  {
    EviThresholds <- seq(0.5,0.9995,by=0.0005)
    
    EFDRgrid <- sapply(EviThresholds, HiddenEFDR, VarThreshold = VarThreshold, Prob = Prob)
    EFNRgrid <- sapply(EviThresholds, HiddenEFNR, VarThreshold = VarThreshold, Prob = Prob)
    
    above <- abs(EFDRgrid - EFDR)
    
    if(sum(!is.na(above)) > 0)
    {
      # Search EFDR closest to the desired value
      EFDRopt <- EFDRgrid[above == min(above, na.rm = TRUE) & !is.na(above)] 
      # If multiple threholds lead to same EFDR, choose the one with the lowest EFDR
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
      OptThreshold <- c(EviThresholds[optimal], EFDRgrid[optimal], EFNRgrid[optimal])
      
      if(OptThreshold[1] < 0.5) 
      {
        message("For the given variance contribution threshold, the evidence threshold 
                that achieves the desired EFDR is below 0.5. By default, the evidence
                threshold will be set at 0.5 (corresponding EFDR/EFNR reported)")
        OptThreshold <- c(EviThresholds[1], EFDRgrid[1], EFNRgrid[1])
      }
      
      # Message when different to desired EFDR is large
      if( abs(OptThreshold[2] - EFDR) > 0.025 )
      {
        message("For the given variance contribution threshold, it is not possible 
                to find an evidence threshold (>0.5) that achieves the desired EFDR level 
                (tolerance +- 0.025). The output of this function reflects the 
                closest possible value. \n")         
      }  
      }
    else
    {
      message("For the given variance contribution threshold, it is not possible 
              to estimate EFDR. By default, the evidence
              threshold will be set at 0.5. \n")    
      OptThreshold <- c(0.5, NA, NA)  
    }
      }
  # If EviThreshold is set a priori
  else
  {
    EFDR = HiddenEFDR(EviThreshold, VarThreshold, Prob)
    EFNR = HiddenEFNR(EviThreshold, VarThreshold, Prob)
    OptThreshold <- c(EviThreshold, EFDR, EFNR)
  }  
  list("EviThresholds" = EviThresholds,
       "EFDRgrid" = EFDRgrid,
       "EFNRgrid" = EFNRgrid,
       "OptThreshold" = OptThreshold)
      }