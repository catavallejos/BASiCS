#' @name BASiCS_VarianceDecomp
#' @aliases BASiCS_VarianceDecomp
#'
#' @title Decomposition of gene expression variability according to BASiCS
#'
#' @description Function to decompose total variability of gene expression into biological and technical components.
#'
#' @param Data an object of class \code{\link[BASiCS]{BASiCS_Data-class}}
#' @param object an object of class \code{\link[BASiCS]{BASiCS_Chain-class}}
#' @param OrderVariable Ordering variable for output. Must take values in \code{c("GeneNames", "BioVarGlobal", "TechVarGlobal", "ShotNoiseGlobal")}.
#' @param Plot If \code{TRUE}, a barplot of the variance decomposition (global and by batches, if any) is generated
#' @param ... Other arguments to be passed to \code{\link[graphics]{barplot}}
#'
#' @return A \code{\link[base]{data.frame}} whose first 4 columns correspond to
#' \describe{
#' \item{\code{GeneName}}{Gene name (as indicated by user)}
#' \item{\code{BioVarGlobal}}{Percentage of variance explained by a biological cell-to-cell heterogeneity component (overall across all cells)}
#' \item{\code{TechVarGlobal}}{Percentage of variance explained by the technical cell-to-cell heterogeneity component (overall across all cells)}
#' \item{\code{ShotNoiseGlobal}}{Percentage of variance explained by the shot noise component (baseline, overall across all cells)}
#' }
#' If more than 1 batch of cells are being analysed, the remaining columns contain the corresponding variance decomposition calculated within each batch.
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
#' @rdname BASiCS_VarianceDecomp

# Change Data class

BASiCS_VarianceDecomp <- function(Data,
                                  object,
                                  OrderVariable = "BioVarGlobal",
                                  Plot = TRUE,
                                  ...)
{
  if(!(OrderVariable %in% c("GeneNames", "BioVarGlobal", "TechVarGlobal", "ShotNoise"))) stop("Invalid 'OrderVariable' value.")
  
  q.bio = ncol(object@delta)
  UniqueBatch = unique(Data@BatchInfo)
  nBatch = length(UniqueBatch)
  
  # Calculating variance decomposition
  VarDecomp = HiddenVarDecomp(Data, object)
  
  # Global values
  BioVarGlobal = apply(VarDecomp$BioVarGlobal, 2, median)
  TechVarGlobal = apply(VarDecomp$TechVarGlobal, 2, median)
  ShotNoiseGlobal = 1-BioVarGlobal-TechVarGlobal
  
  Genes = 1:q.bio
  GeneNames = Data@GeneNames[!Data@Tech]
  
  if(nBatch > 1)
  {
    VarDecompBatch = NULL
    
    for(Batch in 1:nBatch)
    {
      VarDecompBatch = cbind(VarDecompBatch,
                             apply(VarDecomp$BioVarBatch[,,Batch], 2, median),
                             apply(VarDecomp$TechVarBatch[,,Batch], 2, median),
                             1 - apply(VarDecomp$BioVarBatch[,,Batch], 2, median) - apply(VarDecomp$TechVarBatch[,,Batch], 2, median))
    }
    
    colnames(VarDecompBatch) = paste0(rep(c("BioVarBatch", "TechBatch", "ShotNoiseBatch"),nBatch),rep(1:nBatch,each = 3))
    
    out = cbind.data.frame("GeneIndex" = Genes,
                           "GeneNames" = GeneNames,
                           "BioVarGlobal" = BioVarGlobal,
                           "TechVarGlobal" = TechVarGlobal,
                           "ShotNoiseGlobal" = ShotNoiseGlobal,
                           VarDecompBatch,
                           stringsAsFactors = FALSE)
  }
  else
  {
    out = cbind.data.frame("GeneIndex" = Genes,
                           "GeneNames" = GeneNames,
                           "BioVarGlobal" = BioVarGlobal,
                           "TechVarGlobal" = TechVarGlobal,
                           "ShotNoiseGlobal" = ShotNoiseGlobal,
                           stringsAsFactors = FALSE)
  }
  rownames(out) = Genes
  
  # Re-order before output
  if(OrderVariable == "GeneNames") orderVar = GeneNames
  if(OrderVariable == "BioVarGlobal") orderVar = BioVarGlobal
  if(OrderVariable == "TechVarGlobal") orderVar = TechVarGlobal
  if(OrderVariable == "ShotNoiseGlobal") orderVar = 1-BioVarGlobal-TechVarGlobal
  out = out[order(orderVar, decreasing = TRUE),]
  
  if(Plot)
  {
    args <- list(...)
    main = ifelse("main" %in% names(args),args$main, "Overall variance decomposition")
    ylab = ifelse("ylab" %in% names(args),args$ylab, "% of variance")
    beside = ifelse("beside" %in% names(args),args$beside, FALSE)
    if("col" %in% names(args)) {col = args$col} else{col = c("lightblue", "mistyrose", "lightcyan")}
    if("legend" %in% names(args)) {legend = args$legend} else{legend = c("Biological", "Technical", "Shot noise")}
    if("args.legend" %in% names(args)) {args.legend = args$args.legend} else{args.legend = list(x = "bottomright", bg = "white")}
    if("names.arg" %in% names(args)) {names.arg = args$names.arg}
    else
    {
      if(nBatch > 1) {names.arg = c("Overall", paste("Batch ", 1:nBatch))}
      else {names.arg = c("Overall")}
    }
    
    outmat = 100 * matrix(apply(out[,-c(1:2)], 2, mean), nrow = 3, byrow = FALSE)
    barplot(outmat,
            beside = beside, main = main, ylab = ylab,
            col = col, legend = legend,
            args.legend = args.legend,
            names.arg = names.arg)
  }
  
  return(out)
}