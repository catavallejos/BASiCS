#' @name BASiCS_DetectHVG
#' @aliases BASiCS_DetectHVG BASiCS_DetectHVG_LVG
#'
#' @title Detection method for highly and lowly variable genes
#'
#' @description Functions to detect highly and lowly variable genes
#'
#' @param Data an object of class \code{\link[BASiCS]{BASiCS_Data-class}}
#' @param object an object of class \code{\link[BASiCS]{BASiCS_Chain-class}}
#' @param VarThreshold Variance contribution threshold (must be a positive value, between 0 and 1)
#' @param EviThreshold Optional parameter. Evidence threshold (must be a positive value, between 0 and 1)
#' @param OrderVariable Ordering variable for output. Must take values in \code{c("GeneIndex", "Mu", "Delta", "Sigma", "Prob")}.
#' @param Plot If \code{Plot = T} a plot of the gene specific expression level against HVG or LVG is generated.
#' @param ... Graphical parameters (see \code{\link[graphics]{par}}).
#'
#' @return \code{BASiCS_DetectHVG} returns a list of 4 elements:
#' \describe{
#' \item{\code{Table}}{Matrix whose columns contain}
#'    \describe{
#'    \item{\code{Mu}}{Vector of length \code{q.bio}. For each biological gene, posterior median of gene-specific expression levels \eqn{\mu[i]}}
#'    \item{\code{Delta}}{Vector of length \code{q.bio}. For each biological gene, posterior median of gene-specific biological cell-to-cell heterogeneity hyper-parameter \eqn{\delta[i]}}
#'    \item{\code{Sigma}}{Vector of length \code{q.bio}. For each biological gene, proportion of the total variability that is due to a cell-to-cell biological heterogeneity component. }
#'    \item{\code{Prob}}{Vector of length \code{q.bio}. For each biological gene, probability of being highly variable according to the given thresholds.}
#'    \item{\code{HVG}}{Vector of length \code{q.bio}. For each biological gene, indicator of being detected as highly variable according to the given thresholds. }
#'    }
#' \item{\code{EviThreshold}}{Evidence threshold.}
#' \item{\code{EFDR}}{Expected false discovery rate for the given thresholds.}
#' \item{\code{EFNR}}{Expected false negative rate for the given thresholds.}
#' }
#' \code{BASiCS_DetectLVG} produces a similar output, replacing the element \code{HVG} by \code{LVG}, an indicator of a gene being detected as lowly variable according to the given thresholds.
#'
#' @examples
#'
#' # See
#' help(BASiCS_MCMC)
#'
#' @details See vignette
#'
#'
#' @seealso \code{\link[BASiCS]{BASiCS_Chain-class}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#'
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
#'
#' @rdname BASiCS_DetectHVG_LVG
BASiCS_DetectHVG <- function(Data,
                             object,
                             VarThreshold,
                             EviThreshold = NULL,
                             OrderVariable = "Prob",
                             Plot = FALSE,
                             ...)
{
  if(!is(object,"BASiCS_Chain")) stop("'object' is not a BASiCS_Chain class object.")
  if(VarThreshold<0 | VarThreshold>1 | !is.finite(VarThreshold)) stop("Variance contribution thresholds for HVG/LVG detection must be contained in (0,1)")
  if(!is.logical(Plot) | length(Plot)!=1) stop("Please insert T or F for Plot parameter")
  if(!is.null(EviThreshold))
  {
    if(EviThreshold<0 | EviThreshold>1 | !is.finite(EviThreshold))
      stop("Evidence thresholds for HVG and LVG detection must be contained in (0,1) \n For automatic threshold search use EviThreshold = NULL.")
  }
  if(!(OrderVariable %in% c("GeneNames", "Mu", "Delta", "Sigma", "Prob"))) stop("Invalid 'OrderVariable' value")
  Search = F
  if(is.null(EviThreshold)) Search = T
  
  VarDecomp <- HiddenVarDecomp(Data, object)
  Prob <- HiddenProbHVG(VarThreshold = VarThreshold, VarDecomp = VarDecomp)
  
  if(length(EviThreshold) == 0)
  {
    EviThresholds <- seq(0.5,0.9995,by=0.0005)
    
    EFDR <- sapply(EviThresholds, HiddenEFDR, VarThreshold = VarThreshold, Prob = Prob)
    EFNR <- sapply(EviThresholds, HiddenEFNR, VarThreshold = VarThreshold, Prob = Prob)
    
    above<-EFDR>EFNR
    optimal<-which(diff(above)!=0)
    EviThreshold = EviThresholds[optimal]
    if(length(optimal)>0){OptThreshold <- c(EviThreshold, EFDR[optimal], EFNR[optimal])}
    else
    {
      print("It is not possible to find an optimal evidence threshold for the given variance contribution threshold. \n")
      optimal <- round(median(which(abs(EFDR - EFNR) == min(abs(EFDR - EFNR), na.rm = T))))
      if(length(optimal)>0)
      {
        print("Returned value is such that the difference between EFDR and EFNR is minimised.")
        OptThreshold <- c(EviThreshold, EFDR[optimal], EFNR[optimal])
      }
      else
      {
        cat("Numerical issues when computing EFDR and EFNR. Please try a different variance contribution threshold")
        OptThreshold <- rep("Not found",3)
      }
    }
  }
  else
  {
    EFDR = HiddenEFDR(EviThreshold, VarThreshold, Prob)
    EFNR = HiddenEFNR(EviThreshold, VarThreshold, Prob)
    OptThreshold <- c(EviThreshold, EFDR, EFNR)
  }
  
  Sigma <- apply(VarDecomp$BioVarGlobal, 2, median)
  Mu <- apply(object@mu[,1:length(Sigma)], 2, median)
  Delta <- apply(object@delta, 2, median)
  if(OptThreshold[1] == "Not found") {HVG = rep("Not found", length(Sigma))}
  else
  {
    HVG <- ifelse(Prob > OptThreshold[1], TRUE, FALSE);
    HVG <- ifelse(Prob >= 0.5, HVG, FALSE);
  }
  
  qbio = length(Sigma)
  Genes = 1:qbio
  GeneNames = Data@GeneNames[!Data@Tech]
  
  if(Plot)
  {
    args <- list(...)
    
    if(Search)
    {
      par(ask=T)
      
      plot(EviThresholds, EFDR, type = "l", lty = 1, bty = "n", ylab = "Error rate", xlab = "Evidence threshold", ylim = c(0,1))
      lines(EviThresholds, EFNR, lty = 2)
      legend('topleft', c("EFDR", "EFNR"), lty = 1:2, bty = "n")
    }
    
    if("ylim" %in% names(args)) {ylim = args$ylim} else{ylim = c(0, 1)}
    if("xlim" %in% names(args)) {xlim = args$xlim} else{xlim = c(min(Mu),max(Mu))}
    cex = ifelse("cex" %in% names(args),args$cex, 1.5)
    pch = ifelse("pch" %in% names(args),args$pch, 16)
    col = ifelse("col" %in% names(args),args$col, 8)
    bty = ifelse("bty" %in% names(args),args$bty, "n")
    cex.lab = ifelse("cex.lab" %in% names(args),args$cex.lab, 1)
    cex.axis = ifelse("cex.axis" %in% names(args),args$cex.axis, 1)
    cex.main = ifelse("cex.main" %in% names(args),args$cex.main, 1)
    xlab = ifelse("xlab" %in% names(args),args$xlab, expression(mu[i]))
    ylab = ifelse("ylab" %in% names(args),args$ylab, "HVG probability")
    main = ifelse("main" %in% names(args),args$main, "")
    
    plot(Mu, Prob, log="x", pch = pch, ylim = ylim, xlim = xlim, col = col, cex = cex,
         bty = bty, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
         xlab = xlab, ylab = ylab, main = main)
    abline(h = OptThreshold[1], lty = 2, col = "black")
    points(Mu[HVG], Prob[HVG], pch = pch, col = "red", cex = cex)
    
    par(ask=F)
  }
  
  cat(paste(sum(HVG), " genes classified as highly variable using: \n"))
  cat(paste("- Variance contribution threshold = ", round(100*VarThreshold,2), "% \n"))
  cat(paste("- Evidence threshold = ", OptThreshold[1], "\n"))
  cat(paste("- EFDR = ", round(100*OptThreshold[2],2), "% \n"))
  cat(paste("- EFNR = ", round(100*OptThreshold[3],2), "% \n"))
  
  GeneIndex = 1:length(Mu)
  Table = cbind.data.frame("GeneIndex" = Genes,
                           "GeneNames" = GeneNames,
                           "mu" = Mu,
                           "delta" = Delta,
                           "Sigma" = Sigma,
                           "Prob" = Prob,
                           "HVG" = HVG, stringsAsFactors = FALSE)
  rownames(Table) = Genes
  if(OrderVariable == "GeneNames") orderVar = GeneNames
  if(OrderVariable == "Mu") orderVar = Mu
  if(OrderVariable == "Delta") orderVar = Delta
  if(OrderVariable == "Sigma") orderVar = Sigma
  if(OrderVariable == "Prob") orderVar = Prob
  Table = Table[order(orderVar, decreasing = TRUE),]
  
  list("Table" = Table,
       "EviThreshold" = OptThreshold[1], "EFDR" = OptThreshold[2], "EFNR" = OptThreshold[3])
}

#' @name BASiCS_DetectLVG
#' @aliases BASiCS_DetectLVG BASiCS_DetectHVG_LVG
#' @rdname BASiCS_DetectHVG_LVG
BASiCS_DetectLVG <- function(Data,
                             object,
                             VarThreshold,
                             EviThreshold = NULL,
                             OrderVariable = "Prob",
                             Plot = FALSE,
                             ...)
{
  if(!is(object,"BASiCS_Chain")) stop("'object' is not a BASiCS_Chain class object.")
  if(VarThreshold<0 | VarThreshold>1 | !is.finite(VarThreshold)) stop("Variance contribution thresholds for HVG/LVG detection must be contained in (0,1)")
  if(!is.logical(Plot) | length(Plot)!=1) stop("Please insert T or F for Plot parameter")
  if(!is.null(EviThreshold))
  {
    if(EviThreshold<0 | EviThreshold>1 | !is.finite(EviThreshold))
      stop("Evidence thresholds for HVG and LVG detection must be contained in (0,1) \n For automatic threshold search use EviThreshold = NULL.")
  }
  if(!(OrderVariable %in% c("GeneNames", "GeneIndex", "Mu", "Delta", "Sigma", "Prob"))) stop("Invalid 'OrderVariable' value")
  Search = F
  if(is.null(EviThreshold)) Search = T
  
  VarDecomp <- HiddenVarDecomp(Data, object)
  Prob <- HiddenProbLVG(VarThreshold = VarThreshold, VarDecomp = VarDecomp)
  
  if(length(EviThreshold) == 0)
  {
    EviThresholds <- seq(0.5,0.9995,by=0.0005)
    EFDR <- sapply(EviThresholds, HiddenEFDR, VarThreshold = VarThreshold, Prob = Prob)
    EFNR <- sapply(EviThresholds, HiddenEFNR, VarThreshold = VarThreshold, Prob = Prob)
    
    above<-EFDR>EFNR
    optimal<-which(diff(above)!=0)
    EviThreshold = EviThresholds[optimal]
    if(length(optimal)>0){OptThreshold <- c(EviThreshold, EFDR[optimal], EFNR[optimal])}
    else
    {
      print("It is not possible to find an optimal evidence threshold for the given variance contribution threshold. \n")
      optimal <- round(median(which(abs(EFDR - EFNR) == min(abs(EFDR - EFNR), na.rm = T))))
      if(length(optimal)>0)
      {
        print("Returned value is such that the difference between EFDR and EFNR is minimised.")
        OptThreshold <- c(EviThreshold, EFDR[optimal], EFNR[optimal])
      }
      else
      {
        cat("Numerical issues when computing EFDR and EFNR. Please try a different variance contribution threshold")
        OptThreshold <- rep("Not found",3)
      }
    }
  }
  else
  {
    EFDR = HiddenEFDR(EviThreshold, VarThreshold, Prob)
    EFNR = HiddenEFNR(EviThreshold, VarThreshold, Prob)
    OptThreshold <- c(EviThreshold, EFDR, EFNR)
  }
  
  Sigma <- apply(VarDecomp$BioVarGlobal, 2, median)
  Mu <- apply(object@mu[,1:length(Sigma)], 2, median)
  Delta <- apply(object@delta, 2, median)
  if(OptThreshold[1] == "Not found") {LVG = rep("Not found", length(Sigma))}
  else
  {
    LVG <- ifelse(Prob > OptThreshold[1], TRUE, FALSE)
    LVG <- ifelse(Prob >= 0.5, LVG, FALSE);
  }
  
  qbio = length(Sigma)
  Genes = 1:qbio
  GeneNames = Data@GeneNames[!Data@Tech]
  
  if(Plot)
  {
    args <- list(...)
    
    if(Search)
    {
      par(ask=T)
      
      plot(EviThresholds, EFDR, type = "l", lty = 1, bty = "n", ylab = "Error rate", xlab = "Evidence threshold", ylim = c(0,1))
      lines(EviThresholds, EFNR, lty = 2)
      legend('topleft', c("EFDR", "EFNR"), lty = 1:2, bty = "n")
    }
    
    if("ylim" %in% names(args)) {ylim = args$ylim} else{ylim = c(0, 1)}
    if("xlim" %in% names(args)) {xlim = args$xlim} else{xlim = c(min(Mu),max(Mu))}
    cex = ifelse("cex" %in% names(args),args$cex, 1.5)
    pch = ifelse("pch" %in% names(args),args$pch, 16)
    col = ifelse("col" %in% names(args),args$col, 8)
    bty = ifelse("bty" %in% names(args),args$bty, "n")
    cex.lab = ifelse("cex.lab" %in% names(args),args$cex.lab, 1)
    cex.axis = ifelse("cex.axis" %in% names(args),args$cex.axis, 1)
    cex.main = ifelse("cex.main" %in% names(args),args$cex.main, 1)
    xlab = ifelse("xlab" %in% names(args),args$xlab, expression(mu[i]))
    ylab = ifelse("ylab" %in% names(args),args$ylab, "LVG probability")
    main = ifelse("main" %in% names(args),args$main, "")
    
    plot(Mu, Prob, log="x", pch = pch, ylim = ylim, xlim = xlim, col = col, cex = cex,
         bty = bty, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
         xlab = xlab, ylab = ylab, main = main)
    abline(h = OptThreshold[1], lty = 2, col = "black")
    points(Mu[LVG], Prob[LVG], pch = pch, col = "red", cex = cex)
    
    par(ask=F)
    
  }
  
  cat(paste(sum(LVG), " genes classified as lowly variable using: \n"))
  cat(paste("- Variance contribution threshold = ", round(100*VarThreshold,2), "% \n"))
  cat(paste("- Evidence threshold = ", OptThreshold[1], "\n"))
  cat(paste("- EFDR = ", round(100*OptThreshold[2],2), "% \n"))
  cat(paste("- EFNR = ", round(100*OptThreshold[3],2), "% \n"))
  
  GeneIndex = 1:length(Mu)
  Table = cbind.data.frame("GeneIndex" = Genes,
                           "GeneNames" = GeneNames,
                           "mu" = Mu,
                           "delta" = Delta,
                           "Sigma" = Sigma,
                           "Prob" = Prob,
                           "LVG" = LVG, stringsAsFactors = FALSE)
  rownames(Table) = Genes
  
  if(OrderVariable == "GeneNames") orderVar = GeneNames
  if(OrderVariable == "Mu") orderVar = Mu
  if(OrderVariable == "Delta") orderVar = Delta
  if(OrderVariable == "Sigma") orderVar = Sigma
  if(OrderVariable == "Prob") orderVar = Prob
  Table = Table[order(orderVar, decreasing = TRUE),]
  
  list("Table" = Table,
       "EviThreshold" = OptThreshold[1], "EFDR" = OptThreshold[2], "EFNR" = OptThreshold[3])
}
