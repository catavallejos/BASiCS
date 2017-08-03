#### Helper functions for differential testing

HiddenTailProbUpDV<-function(chain,threshold){return(sum(chain>threshold)/length(chain))}
HiddenTailProbLowDV<-function(chain,threshold){return(sum(chain<threshold)/length(chain))}

HiddenProbDE<-function(
  chain, # MCMC chain for log-fold change parameter (in absolute value)
  tol) # Minimum tolerance difference for log-fold change
{  
  return(apply(chain, 2, HiddenTailProbUpDV, threshold = tol))
}

HiddenEFDRDV<-function(
  EviThreshold, # Evidence threshold (it must be contained in (0,1))
  Prob) # Probability of changes in expression (for a given tolerance)
{  
  return(sum((1-Prob)*I(Prob > EviThreshold))/sum(I(Prob > EviThreshold)))
}

HiddenEFNRDV<-function(
  EviThreshold, # Evidence threshold (it must be contained in (0,1))
  Prob) # Probability of changes in expression (for a given tolerance) 
{
  return(sum(Prob*I(Prob <= EviThreshold))/sum(I(Prob <= EviThreshold)))
}

#### Helper functions for detecting highly and lowly variable genes

HiddenTailProbUp<-function(chain,threshold){return(sum(chain>threshold)/length(chain))}
HiddenTailProbLow<-function(chain,threshold){return(sum(chain<threshold)/length(chain))}

HiddenProbHVG<-function(
  VarThreshold, # Variance contribution threshold for HVG/LVG detection. Value must be between 0 and 1.
  VarDecomp) # Output of the variance decomposition obtained using BASiCS_Variance
{
  return(apply(VarDecomp$BioVarGlobal, 2, HiddenTailProbUp, threshold = VarThreshold))
}

HiddenProbLVG<-function(
  VarThreshold, # Variance contribution threshold for HVG/LVG detection. Value must be between 0 and 1.
  VarDecomp) # Output of the variance decomposition obtained using BASiCS_Variance
{
  return(apply(VarDecomp$BioVarGlobal, 2, HiddenTailProbLow, threshold = VarThreshold))
}

HiddenEFDR<-function(
  EviThreshold, # Evidence threshold (it must be contained in (0,1))
  VarThreshold, # Variance contribution threshold choosen by the user (it must be contained in (0,1))
  Prob) # Output of the variance decomposition obtained using BASiCS_Variance
{
  return(sum((1-Prob)*I(Prob>EviThreshold))/sum(I(Prob>EviThreshold)))
}

HiddenEFNR<-function(
  EviThreshold, # Evidence threshold (it must be contained in (0,1))
  VarThreshold, # Variance contribution threshold choosen by the user (it must be contained in (0,1))
  Prob) # Output of the variance decomposition obtained using BASiCS_Variance
{
  return(sum(Prob*I(EviThreshold>=Prob))/sum(I(EviThreshold>=Prob)))
}

HiddenHeaderDetectHVG_LVG <- function(Data,
                                      object,
                                      VarThreshold,
                                      EviThreshold = NULL,
                                      EFDR = 0.05, 
                                      OrderVariable = "Prob",
                                      Plot = FALSE)
{
  if(!is(Data,"SummarizedExperiment")) stop("'Data' is not a SummarizedExperiment class object. Please use the 'newBASiCS_Data' function to create a SummarizedExperiment object.")
  if(!is(object,"BASiCS_Chain")) stop("'object' is not a BASiCS_Chain class object.")
  if(VarThreshold<0 | VarThreshold>1 | !is.finite(VarThreshold)) stop("Variance contribution thresholds for HVG/LVG detection must be contained in (0,1)")
  if(!is.logical(Plot) | length(Plot)!=1) stop("Please insert TRUE or FALSE for Plot parameter")
  if(!is.null(EviThreshold))
  {
    if(EviThreshold<0 | EviThreshold>1 | !is.finite(EviThreshold))
      stop("Evidence thresholds for HVG and LVG detection must be contained in (0,1) \n For automatic threshold search use EviThreshold = NULL.")
  }
  if(!(OrderVariable %in% c("GeneNames", "Mu", "Delta", "Sigma", "Prob"))) stop("Invalid 'OrderVariable' value")
  if(is.null(EviThreshold)) {message(paste("Posterior probability threshold to be defined by EFDR = ", 100*EFDR, "% (+-2.5% tolerance) ..."))}
}

HiddenPlot1DetectHVG_LVG <- function(EviThresholds, 
                                     EFDRgrid,
                                     EFNRgrid,
                                     OptThreshold,
                                     EviThreshold,
                                     EFDR)
{
  plot(EviThresholds, EFDRgrid, type = "l", lty = 1, bty = "n", lwd = 2,  
       ylab = "Error rate", xlab = "Evidence threshold", ylim = c(0,1))
  lines(EviThresholds, EFNRgrid, lty = 2, lwd = 2)
  if(length(EviThreshold) == 0) {abline(h = EFDR, col = "blue", lwd = 2, lty = 1)}
  abline(v = OptThreshold[1], col = "red", lwd = 2, lty = 1)
  if(length(EviThreshold) == 0)
  {
    legend('topleft', c("EFDR", "EFNR", "Target EFDR"), lty = c(1:2,1), 
           col = c("black", "black", "blue"), bty = "n", lwd = 2)  
  }
}

HiddenPlot2DetectHVG_LVG <- function(args, 
                                     Task, 
                                     Mu, 
                                     Prob,
                                     OptThreshold,
                                     Hits)
{
  if("ylim" %in% names(args)) {ylim = args$ylim} else{ylim = c(0, 1)}
  if("xlim" %in% names(args)) {xlim = args$xlim} else{xlim = c(min(Mu),max(Mu))}
  cex = ifelse("cex" %in% names(args),args$cex, 1.5)
  pch = ifelse("pch" %in% names(args),args$pch, 16)
  col = ifelse("col" %in% names(args),args$col, 8)
  bty = ifelse("bty" %in% names(args),args$bty, "n")
  cex.lab = ifelse("cex.lab" %in% names(args),args$cex.lab, 1)
  cex.axis = ifelse("cex.axis" %in% names(args),args$cex.axis, 1)
  cex.main = ifelse("cex.main" %in% names(args),args$cex.main, 1)
  xlab = ifelse("xlab" %in% names(args),args$xlab, "Mean expression")
  if(Task == "HVG") {ylab = ifelse("ylab" %in% names(args),args$ylab, "HVG probability")}
  else {ylab = ifelse("ylab" %in% names(args),args$ylab, "LVG probability")}
  main = ifelse("main" %in% names(args),args$main, "")
  
  plot(Mu, Prob, log="x", pch = pch, ylim = ylim, xlim = xlim, col = col, cex = cex,
       bty = bty, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
       xlab = xlab, ylab = ylab, main = main)
  abline(h = OptThreshold[1], lty = 2, col = "black")
  points(Mu[Hits], Prob[Hits], pch = pch, col = "red", cex = cex)
}



