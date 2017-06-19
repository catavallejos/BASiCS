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

#' @name BASiCS_D_TestDE
#' 
#' @title Detection of genes with changes in expression
#' 
#' @description Function to assess changes in expression (mean and over-dispersion)
#' 
#' @param Data an object of class \code{\link[BASiCS]{BASiCS_D_Data-class}}
#' @param object an object of class \code{\link[BASiCS]{BASiCS_D_Chain-class}}
#' @param GeneNames Vector containing gene names to be used in results table (argument to be removed as 'GeneNames' will be an slot of `BASiCS_D_Data` object)
#' @param EpsilonM Minimum fold change tolerance threshold for detecting changes in overall expression (must be a positive real number)
#' @param EpsilonD Minimum fold change tolerance threshold for detecting changes in cell-to-cell biological over dispersion (must be a positive real number)
#' @param EviThresholdM Optional parameter. Evidence threshold for detecting changes in overall expression (must be a positive value, between 0 and 1)
#' @param EviThresholdD Optional parameter. Evidence threshold for detecting changes in cell-to-cell biological over dispersion (must be a positive value, between 0 and 1)
#' @param OrderVariable Ordering variable for output. Must take values in \code{c("GeneIndex", "GeneNames", "ProbDiffExp", "ProbDiffOverDisp")}.
#' @param GroupLabelRef Label assigned to reference group. Default: \code{GroupLabelRef = "Ref"}
#' @param GroupLabelTest Label assigned to reference group. Default: \code{GroupLabelRef = "Test"}
#' @param Plot If \code{Plot = T}, error rates control rates and volcano plots are generated.  
#' @param OffSet Optional argument to remove a fix offset effect (if not previously removed from the MCMC chains). This argument will be removed shorly, once offset removal is built as an internal step. 
#' @param EFDR_M Target for expected false discovery rate related to the comparison of means (default = 0.05)
#' @param EFDR_D Target for expected false discovery rate related to the comparison of dispersions (default = 0.05)
#' @param GenesSelect Optional argument to provide a user-defined list of genes to be considered for the comparison (default = NULL). When used, this argument must be a vector of 'TRUE' (include gene) / 'FALSE' (exclude gene) indicator, with the same length as the number of intrinsic genes and following the same order as how genes are displayed in the table of counts.  This argument is necessary in order to have a meaningful EFDR calibration when the user decides to exclude some genes from the comparison. 
#' @param ... Graphical parameters (see \code{\link[graphics]{par}}).
#' 
#' @return \code{BASiCS_D_TestDE} returns a list of 4 elements:
#' \describe{
#' \item{\code{Table}}{A \code{\link[base]{data.frame}} containing the results of the differential expression test}
#'    \describe{
#'    \item{\code{Mu}}{Vector of length \code{q.bio}. For each biological gene, posterior median of gene-specific expression levels \eqn{\mu[i]}}
#'    \item{\code{Delta}}{Vector of length \code{q.bio}. For each biological gene, posterior median of gene-specific biological cell-to-cell heterogeneity hyper-parameter \eqn{\delta[i]}}
#'    \item{\code{Sigma}}{Vector of length \code{q.bio}. For each biological gene, proportion of the total variability that is due to a cell-to-cell biological heterogeneity component. }
#'    \item{\code{Prob}}{Vector of length \code{q.bio}. For each biological gene, probability of being highly variable according to the given thresholds.}
#'    \item{\code{HVG}}{Vector of length \code{q.bio}. For each biological gene, indicator of being detected as highly variable according to the given thresholds. }
#'    }
#' \item{\code{EviThreshold}}{Evidence thresholds.}
#' \item{\code{EFDR}}{Expected false discovery rate for the given thresholds.}
#' \item{\code{EFNR}}{Expected false negative rate for the given thresholds.}
#' }
#'  
#' @examples
#' 
#' # See vignette
#' 
#' @details See vignette
#' 
#' 
#' @seealso \code{\link[BASiCS]{BASiCS_D_Chain-class}} 
#' 
#' @author Catalina A. Vallejos \email{cnvallej@uc.cl}
#' 
#' @rdname BASiCS_D_TestDE

# Set default OffSet to TRUE + warning when Offset = FALSE
# Leave the .cl email
# Change input to 2 SummarizedExperiment and 2 MCMC
# Remove Data from input

BASiCS_D_TestDE <- function(Data, 
                            object,
                            GeneNames,
                            EpsilonM = 0.10,
                            EpsilonD = 0.10,
                            EviThresholdM = NULL,
                            EviThresholdD = NULL,
                            OrderVariable = "ProbDiffExp",
                            GroupLabelRef = "Ref",
                            GroupLabelTest = "Test",
                            Plot = FALSE, 
                            OffSet = FALSE, 
                            EFDR_M = 0.05,
                            EFDR_D = 0.05,
                            GenesSelect = NULL, 
                            ...)
{
  # Checking validity of input arguments
  if(!is(Data,"BASiCS_D_Data")) stop("'Data' is not a 'BASiCS_D_Data' class object.")    
  if(!is(object,"BASiCS_D_Chain")) stop("'object' is not a 'BASiCS_D_Chain' class object.")    
  if(EpsilonM < 0 | !is.finite(EpsilonM)) stop("Minimum tolerance of fold-change 'EpsilonM' must be a positive real value")
  if(EpsilonD < 0 | !is.finite(EpsilonD)) stop("Minimum tolerance of fold-change 'EpsilonD' must be a positive real value")
  if(!is.logical(Plot) | length(Plot)!=1) stop("Please insert T or F for Plot parameter")
  if(!is.null(EviThresholdM) | !is.null(EviThresholdD))
  {
    if(EviThresholdM < 0 | EviThresholdM > 1 | !is.finite(EviThresholdM)) 
      stop("Evidence threshold 'EviThresholdM' must be contained in (0,1) \n For automatic threshold search use EviThresholdM = NULL.")    
    if(EviThresholdD < 0 | EviThresholdD > 1 | !is.finite(EviThresholdD)) 
      stop("Evidence threshold 'EviThresholdD' must be contained in (0,1) \n For automatic threshold search use EviThresholdD = NULL.")    
  }
  if(!(OrderVariable %in% c("GeneIndex", "GeneNames", "ProbDiffExp", "ProbDiffOverDisp"))) stop("Invalid 'OrderVariable' value")
  if(!is.character(GroupLabelRef) | length(GroupLabelRef) > 1) stop("Invalid value for 'GroupLabelRef'")
  if(!is.character(GroupLabelTest) | length(GroupLabelTest) > 1) stop("Invalid value for 'GroupLabelTest'")
  if(!is.null(GenesSelect) & (length(GenesSelect) != length(GeneNames))) stop("Invalid value for 'GenesSelect'")
  if(!is.null(GenesSelect) & !is.logical(GenesSelect)) stop("Invalid value for 'GenesSelect'")
  # if(!is.null(GenesSelect) & (length(GenesSelect) == sum(!GenesSelect))) stop("Invalid value for 'GenesSelect'")
  
  nTest = ncol(Data@CountsTest)
  nRef = ncol(Data@CountsRef)
  n = nTest + nRef
  # Changes in overall expression
  if(OffSet)
  {
    # With offset correction
    ChainMuRefOffSet = object@muRef / rowSums(object@muRef)
    ChainMuTestOffSet = object@muTest / rowSums(object@muTest)
    MedianMuRefOffSet = apply(ChainMuRefOffSet, 2, median)
    MedianMuTestOffSet = apply(ChainMuTestOffSet, 2, median)
    
    ChainTau = log(ChainMuTestOffSet / ChainMuRefOffSet )
    MedianTau = apply(ChainTau, 2, median)
  }
  else
  {
    ChainTau = log(object@muTest / object@muRef)
    MedianTau = apply(ChainTau, 2, median)  
  }
  
  if(EpsilonM > 0) {ProbM = HiddenProbDE(chain = abs(ChainTau), tol = EpsilonM)}
  else 
  {
    ProbM_aux = HiddenProbDE(chain = ChainTau, tol = 0)
    ProbM = 2*pmax(ProbM_aux, 1-ProbM_aux) - 1
  }
  if(is.null(EviThresholdM))
  {
    EviThresholds <- seq(0.5,0.9995,by=0.00025)
    
    if(is.null(GenesSelect))
    {
      EFDR <- sapply(EviThresholds, HiddenEFDRDV, Prob = ProbM)
      EFNR <- sapply(EviThresholds, HiddenEFNRDV, Prob = ProbM)      
    }
    else
    {
      EFDR <- sapply(EviThresholds, HiddenEFDRDV, Prob = ProbM[GenesSelect])
      EFNR <- sapply(EviThresholds, HiddenEFNRDV, Prob = ProbM[GenesSelect])        
    }
    
    optimal = round(median(which(abs(EFDR - EFDR_M) == min(abs(EFDR - EFDR_M)))))
    OptThresholdM <- c(EviThresholds[optimal], EFDR[optimal], EFNR[optimal])
    
    EviThresholdM = OptThresholdM[1]
    if(is.na(EviThresholdM)) 
    { 
      cat("EFDR calibration failed. Probability threshold automatically set equal to 0.90 \n")
      EviThresholdM = 0.90
    }
  }
  else
  {
    if(is.null(GenesSelect))
    {
      EFDR = HiddenEFDRDV(EviThresholdM, ProbM)
      EFNR = HiddenEFNRDV(EviThresholdM, ProbM)
    }
    else
    {
      EFDR = HiddenEFDRDV(EviThresholdM, ProbM[GenesSelect])
      EFNR = HiddenEFNRDV(EviThresholdM, ProbM[GenesSelect])      
    }
    OptThresholdM <- c(EviThresholdM, EFDR, EFNR)
  }  
  
  # Changes in cell-to-cell biological over dispersion
  ChainOmega = log(object@deltaTest / object@deltaRef)
  MedianOmega = apply(ChainOmega, 2, median)
  if(EpsilonD > 0) {ProbD = HiddenProbDE(chain = abs(ChainOmega), tol = EpsilonD)}
  else 
  {
    ProbD_aux = HiddenProbDE(chain = ChainOmega, tol = 0)
    ProbD = 2*pmax(ProbD_aux, 1-ProbD_aux) - 1
  }
  if(is.null(EviThresholdD))
  {
    EviThresholds <- seq(0.5,0.9995,by=0.00025)
    
    if(is.null(GenesSelect))
    {
      EFDR <- sapply(EviThresholds, HiddenEFDRDV, Prob = ProbD)
      EFNR <- sapply(EviThresholds, HiddenEFNRDV, Prob = ProbD)
    }
    else
    {
      EFDR <- sapply(EviThresholds, HiddenEFDRDV, Prob = ProbD[GenesSelect])
      EFNR <- sapply(EviThresholds, HiddenEFNRDV, Prob = ProbD[GenesSelect])
    }
    
    optimal = round(median(which(abs(EFDR - EFDR_D) == min(abs(EFDR - EFDR_D)))))
    OptThresholdD <- c(EviThresholds[optimal], EFDR[optimal], EFNR[optimal])
    
    EviThresholdD = OptThresholdD[1]
    if(is.na(EviThresholdD)) 
    { 
      cat("EFDR calibration failed. Probability threshold automatically set equal to 0.90 \n")
      EviThresholdD = 0.90
    }
  }
  else
  {
    if(is.null(GenesSelect))
    {
      EFDR = HiddenEFDRDV(EviThresholdD, ProbD)
      EFNR = HiddenEFNRDV(EviThresholdD, ProbD)
    }
    else
    {
      EFDR = HiddenEFDRDV(EviThresholdD, ProbD[GenesSelect])
      EFNR = HiddenEFNRDV(EviThresholdD, ProbD[GenesSelect])      
    }
    OptThresholdD <- c(EviThresholdD, EFDR, EFNR)
  } 
  
  # Differential expression results
  TestPlusM = which(ProbM > EviThresholdM & MedianTau > 0)
  RefPlusM = which(ProbM > EviThresholdM & MedianTau < 0)  
  ResultDiffExp = rep("NoDiff", length(MedianTau))
  ResultDiffExp[TestPlusM] = paste0(GroupLabelTest,"+")
  ResultDiffExp[RefPlusM] = paste0(GroupLabelRef,"+") 
  if(!is.null(GenesSelect)) {ResultDiffExp[!GenesSelect] = "ExcludedByUser"}
  
  # Differential cell-to-cell biological over dispersion results
  TestPlusD = which(ProbD > EviThresholdD & MedianOmega > 0)
  RefPlusD = which(ProbD > EviThresholdD & MedianOmega < 0)
  ResultDiffOverDisp = rep("NoDiff", length(MedianTau))
  ResultDiffOverDisp[TestPlusD] = paste0(GroupLabelTest,"+")
  ResultDiffOverDisp[RefPlusD] = paste0(GroupLabelRef,"+")
  if(!is.null(GenesSelect)) {ResultDiffOverDisp[!GenesSelect] = "ExcludedByUser"}
  
  # Gene-specific parameters' summaries
  MedianMuTest = apply(object@muTest, 2, median)
  MedianMuRef = apply(object@muRef, 2, median)  
  MedianDeltaTest = apply(object@deltaTest, 2, median)
  MedianDeltaRef = apply(object@deltaRef, 2, median)
  MuBase=(MedianMuRef * nRef + MedianMuTest * nTest)/n
  DeltaBase=(MedianDeltaRef * nRef + MedianDeltaTest * nTest)/n
  
  FullList = cbind.data.frame("GeneNames" = GeneNames,
                              "ExpOverall"= round(as.numeric(MuBase),3),
                              "ExpTest"= round(as.numeric(MedianMuTest),3),
                              "ExpRef"= round(as.numeric(MedianMuRef),3),
                              "ExpFC"= round(as.numeric(exp(MedianTau)),3),
                              "ExpLogFC"= round(as.numeric(MedianTau),3),
                              "ProbDiffExp"= round(as.numeric(ProbM),3),
                              "ResultDiffExp" = ResultDiffExp,
                              "OverDispOverall"= round(as.numeric(DeltaBase),3),
                              "OverDispTest"= round(as.numeric(MedianDeltaTest),3),
                              "OverDispRef"= round(as.numeric(MedianDeltaRef),3),  
                              "OverDispFC"= round(as.numeric(exp(MedianOmega)),3), 
                              "OverDispLogFC"= round(as.numeric(MedianOmega),3),
                              "ProbDiffOverDisp"= round(as.numeric(ProbD),3),
                              "ResultDiffOverDisp" = ResultDiffOverDisp,
                              stringsAsFactors = FALSE)
  
  GeneIndex = 1:length(MuBase)
  
  if(OrderVariable == "GeneIndex") orderVar = GeneIndex
  if(OrderVariable == "GeneNames") orderVar = GeneNames
  if(OrderVariable == "ProbDiffExp") orderVar = ProbM
  if(OrderVariable == "ProbDiffOverDisp") orderVar = ProbD
  
  FullList = FullList[order(orderVar, decreasing = TRUE),]
  
  if(!is.null(GenesSelect))
  {
    cat("--------------------------------------------------------------------- \n")
    cat(paste("The user excluded ", sum(!GenesSelect), " genes from the comparison. \n"))
    cat("These genes are marked as 'ExcludedByUser' in the results table and excluded from EFDR calibration. \n")
    cat("--------------------------------------------------------------------- \n")    
  }
  
  nPlusTestM = sum(ResultDiffExp == paste0(GroupLabelTest,"+"))
  nPlusRefM = sum(ResultDiffExp == paste0(GroupLabelRef,"+"))
  nPlusTestD = sum(ResultDiffOverDisp == paste0(GroupLabelTest,"+"))
  nPlusRefD = sum(ResultDiffOverDisp == paste0(GroupLabelRef,"+"))
  
  cat("--------------------------------------------------------------------- \n")
  cat(paste(nPlusTestM + nPlusRefM, " genes with a change on the overall expression:  \n"))
  cat(paste("- Higher expression in ",GroupLabelTest,"group:", nPlusTestM, "\n"))
  cat(paste("- Higher expression in ",GroupLabelRef,"group:", nPlusRefM, "\n"))
  cat(paste("- Fold change tolerance = ", round(100*EpsilonM,2), "% \n"))    
  cat(paste("- Evidence threshold = ", OptThresholdM[1], "\n"))
  cat(paste("- EFDR = ", round(100*OptThresholdM[2],2), "% \n"))  
  cat(paste("- EFNR = ", round(100*OptThresholdM[3],2), "% \n"))  
  cat("--------------------------------------------------------------------- \n")
  cat("\n")
  cat("--------------------------------------------------------------------- \n")
  cat(paste(nPlusTestD + nPlusRefD, " genes with a change on the cell-to-cell biological over dispersion:  \n"))
  cat(paste("- Higher over dispersion in ",GroupLabelTest,"group:", nPlusTestD, "\n"))
  cat(paste("- Higher over dispersion in ",GroupLabelRef,"group:", nPlusRefD, "\n"))
  cat(paste("- Fold change tolerance = ", round(100*EpsilonD,2), "% \n"))    
  cat(paste("- Evidence threshold = ", OptThresholdD[1], "\n"))
  cat(paste("- EFDR = ", round(100*OptThresholdD[2],2), "% \n"))  
  cat(paste("- EFNR = ", round(100*OptThresholdD[3],2), "% \n"))    
  cat("--------------------------------------------------------------------- \n")
  
  if(Plot)
  {    
    args <- list(...)
    
    #    if(Search)
    #    {      
    #      par(ask=T)
    #      
    #      plot(EviThresholds, EFDR, type = "l", lty = 1, bty = "n", ylab = "Error rate", xlab = "Evidence threshold", ylim = c(0,1))
    #      lines(EviThresholds, EFNR, lty = 2)      
    #      legend('topleft', c("EFDR", "EFNR"), lty = 1:2, bty = "n")
    #    }
    
    #    if("ylim" %in% names(args)) {ylim = args$ylim} else{ylim = c(0, 1)}
    #    if("xlim" %in% names(args)) {xlim = args$xlim} else{xlim = c(min(Mu),max(Mu))}
    #    cex = ifelse("cex" %in% names(args),args$cex, 1.5)
    #    pch = ifelse("pch" %in% names(args),args$pch, 16)
    #    col = ifelse("col" %in% names(args),args$col, 8)
    #    bty = ifelse("bty" %in% names(args),args$bty, "n")
    #    cex.lab = ifelse("cex.lab" %in% names(args),args$cex.lab, 1)
    #    cex.axis = ifelse("cex.axis" %in% names(args),args$cex.axis, 1)
    #    cex.main = ifelse("cex.main" %in% names(args),args$cex.main, 1) 
    #    xlab = ifelse("xlab" %in% names(args),args$xlab, expression(mu[i]))
    #    ylab = ifelse("ylab" %in% names(args),args$ylab, "HVG probability")
    #    main = ifelse("main" %in% names(args),args$main, "") 
    
    #    plot(Mu, Prob, log="x", pch = pch, ylim = ylim, xlim = xlim, col = col, cex = cex,
    #         bty = bty, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
    #         xlab = xlab, ylab = ylab, main = main)
    #    abline(h = OptThreshold[1], lty = 2, col = "black")
    #    points(Mu[HVG], Prob[HVG], pch = pch, col = "red", cex = cex)
    
    #    par(ask=F)
  }
  
  list("Table" = FullList, 
       "DiffExpSummary" = list("EviThreshold" = OptThresholdM[1],
                               "EFDR" = OptThresholdM[2], 
                               "EFNR" = OptThresholdM[3]),
       "DiffOverDispSummary" = list("EviThreshold" = OptThresholdD[1],
                                    "EFDR" = OptThresholdD[2], 
                                    "EFNR" = OptThresholdD[3]))
}