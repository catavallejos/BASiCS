#' @name BASiCS_VarianceDecomp
#' @aliases BASiCS_VarianceDecomp
#'
#' @title Decomposition of gene expression variability according to BASiCS
#'
#' @description Function to decompose total variability of gene 
#' expression into biological and technical components.
#'
#' @param Chain an object of class \code{\linkS4class{BASiCS_Chain}} 
#' @param OrderVariable Ordering variable for output. 
#' Possible values: \code{'GeneName'}, \code{'BioVarGlobal'},
#'  \code{'TechVarGlobal'} and \code{'ShotNoiseGlobal'}. 
#' Default: \code{OrderVariable = "BioVarGlobal"}.
#' @param Plot If \code{TRUE}, a barplot of the variance decomposition 
#' (global and by batches, if any) is generated. Default: \code{Plot = TRUE}.
#' @param ... Other arguments to be passed to \code{\link[graphics]{barplot}}
#'
#' @return A \code{\link[base]{data.frame}} whose first 4 columns correspond to
#' \describe{
#' \item{\code{GeneName}}{Gene name (as indicated by user)}
#' \item{\code{BioVarGlobal}}{Percentage of variance explained by a biological 
#'                            component (overall across all cells)}
#' \item{\code{TechVarGlobal}}{Percentage of variance explained by the technical 
#'                             component (overall across all cells)}
#' \item{\code{ShotNoiseGlobal}}{Percentage of variance explained by the shot 
#'                               noise component (baseline Poisson noise, 
#'                               overall across all cells)}
#' }
#' If more than 1 batch of cells are being analysed, the remaining columns 
#' contain the corresponding variance decomposition calculated within each batch.
#'
#' @examples
#'
#' # For illustration purposes we load a built-in 'BASiCS_Chain' object 
#' # (obtained using the 'BASiCS_MCMC' function)
#' data(ChainSC)
#' 
#' VD <- BASiCS_VarianceDecomp(ChainSC)
#'
#' @details See vignette
#'
#' @seealso \code{\linkS4class{BASiCS_Chain}} 
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#'
#' @references 
#' 
#' Vallejos, Marioni and Richardson (2015). PLoS Computational Biology. 
#'
#' @rdname BASiCS_VarianceDecomp
#' @export
BASiCS_VarianceDecomp <- function(Chain, 
                                  OrderVariable = "BioVarGlobal", 
                                  Plot = TRUE, ...) 
{
  if (!(OrderVariable %in% c("GeneName", "BioVarGlobal", 
                             "TechVarGlobal", "ShotNoise"))) 
    stop("Invalid 'OrderVariable' value.")
  if (!is(Chain, "BASiCS_Chain")) 
    stop("'Chain' is not a BASiCS_Chain class object.")
    
  q.bio <- ncol(Chain@parameters$delta)
  UniqueBatch <- colnames(Chain@parameters$theta)
  nBatch <- length(UniqueBatch)
    
  # Calculating variance decomposition
  VarDecomp <- HiddenVarDecomp(Chain)
    
  # Global values
  BioVarGlobal <- matrixStats::colMedians(VarDecomp$BioVarGlobal)
  TechVarGlobal <- matrixStats::colMedians(VarDecomp$TechVarGlobal)
  ShotNoiseGlobal <- 1 - BioVarGlobal - TechVarGlobal
    
  Genes <- seq_len(q.bio)
  GeneName <- colnames(Chain@parameters$mu)
    
  if (nBatch > 1) 
  {
    VarDecompBatch <- NULL
        
    for (Batch in seq_len(nBatch)) 
    {
      BioVarAux <- matrixStats::colMedians(VarDecomp$BioVarBatch[, , Batch])
      TechVarAux <- matrixStats::colMedians(VarDecomp$TechVarBatch[, , Batch])
      VarDecompBatch <- cbind(VarDecompBatch, BioVarAux, 
                              TechVarAux, 1 - BioVarAux - TechVarAux)
    }
        
    colnames(VarDecompBatch) <- paste0(rep(c("BioVarBatch", 
                                             "TechBatch", 
                                             "ShotNoiseBatch"), 
                                              nBatch), 
                                       rep(seq_len(nBatch), each = 3))
        
    out <- cbind.data.frame(GeneIndex = Genes, GeneName = GeneName, 
                            BioVarGlobal = BioVarGlobal, 
                            TechVarGlobal = TechVarGlobal, 
                            ShotNoiseGlobal = ShotNoiseGlobal, 
                            VarDecompBatch, stringsAsFactors = FALSE)
  } 
  else 
  {
    out <- cbind.data.frame(GeneIndex = Genes, GeneName = GeneName, 
                            BioVarGlobal = BioVarGlobal, 
                            TechVarGlobal = TechVarGlobal, 
                            ShotNoiseGlobal = ShotNoiseGlobal, 
                            stringsAsFactors = FALSE)
  }
  rownames(out) <- Genes
    
  # Re-order before output
  if (OrderVariable == "GeneName") { orderVar <- GeneName }
  if (OrderVariable == "BioVarGlobal") { orderVar <- BioVarGlobal }
  if (OrderVariable == "TechVarGlobal") { orderVar <- TechVarGlobal }
  if (OrderVariable == "ShotNoiseGlobal") 
    orderVar <- 1 - BioVarGlobal - TechVarGlobal
  
  out <- out[order(orderVar, decreasing = TRUE), ]
    
  if (Plot) 
  {
    args <- list(...)
    main <- ifelse("main" %in% names(args), 
                   args$main, "Overall variance decomposition")
    ylab <- ifelse("ylab" %in% names(args), args$ylab, "% of variance")
    beside <- ifelse("beside" %in% names(args), args$beside, FALSE)
    
    if ("col" %in% names(args)) { col <- args$col } 
    else { col <- c("lightblue", "mistyrose", "lightcyan") }
    
    if ("legend" %in% names(args)) { legend = args$legend } 
    else { legend <- c("Biological", "Technical", "Shot noise") }
    
    if ("args.legend" %in% names(args)) { args.legend = args$args.legend } 
    else { args.legend <- list(x = "bottomright", bg = "white") }
    
    if ("names.arg" %in% names(args)) { names.arg <- args$names.arg } 
    else 
    {
      if (nBatch > 1) 
      { 
        names.arg <- c("Overall", paste("Batch ", seq_len(nBatch))) 
      } 
      else { names.arg <- c("Overall") } 
    }
        
    outmat <- 100 * matrix(apply(out[, -c(1,2)], 2, mean), 
                           nrow = 3, byrow = FALSE)
    barplot(outmat, beside = beside, main = main, 
            ylab = ylab, col = col, 
            legend = legend, args.legend = args.legend, 
            names.arg = names.arg)
  }
  return(out)
}
