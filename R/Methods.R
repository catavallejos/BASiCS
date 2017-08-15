# Remove functions that are not necessary anymore

##########################################################################
# New generic functions
##########################################################################
setGeneric("displayChainBASiCS", function(object, ...) standardGeneric("displayChainBASiCS"))
setGeneric("displaySummaryBASiCS", function(object, ...) standardGeneric("displaySummaryBASiCS"))
#setGeneric("displayOffsetBASiCS", function(object, ...){})

##########################################################################
# Methods for BASiCS_Chain objects
##########################################################################

#' @name BASiCS_Chain-methods
#' @aliases show,BASiCS_Chain-method
#' 
#' @title 'show' method for BASiCS_Chain objects
#' 
#' @description 'show' method for \code{\link[BASiCS]{BASiCS_Chain-class}} objects.
#' 
#' @param object A \code{BASiCS_Chain} object.
#' 
#' @return Prints a summary of the properties of \code{object}.
#' 
#' @examples
#' 
#' Data = makeExampleBASiCS_Data()
#' Chain <- BASiCS_MCMC(Data, N = 50, Thin = 2, Burn = 2)
#' 
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' 
#' @references 
#' 
#' Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. PLoS Computational Biology.
#' 
#' @rdname BASiCS_Chain-methods
setMethod("show",
          signature = "BASiCS_Chain",
          definition = function(object){
            N = nrow(object@mu)
            q.bio = ncol(object@delta)
            n = ncol(object@phi)
            nBatch = ncol(object@theta)
            message("An object of class ", class(object), "\n", sep = "")
            message(" ", N," MCMC samples.\n", sep = "")
            if(nBatch > 1)
            {
              message(" Dataset contains ", q.bio, " biological genes and ", n, " cells (",nBatch," batches). \n", sep="")
            }
            else{message(" Dataset contains ", q.bio, " biological genes and ", n, " cells (1 batch). \n", sep="")}
            message(" Elements (slots): mu, delta, phi, s, nu and theta.\n")
          })

#' @name Summary
#' @aliases Summary Summary,BASiCS_Chain-method
#' 
#' @docType methods
#' @rdname Summary-BASiCS_Chain-method
#' 
#' @title 'Summary' method for BASiCS_Chain objects
#' 
#' @description For each of the BASiCS parameters (see Vallejos et al 2015), 
#' \code{Summary} returns the corresponding postior medians and limits of the high posterior
#' density interval with probabilty equal to \code{prob}.
#' 
#' @param x A \code{BASiCS_Chain} object.
#' @param prob \code{prob} argument for \code{\link[coda]{HPDinterval}} function. 
#' 
#' @return An object of class \code{\link[BASiCS]{BASiCS_Summary-class}}. 
#' 
#' @examples 
#' 
#' help(BASiCS_MCMC)
#' 
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' 
#' @references 
#' 
#' Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. PLoS Computational Biology.
setMethod("Summary",
          signature = "BASiCS_Chain",
          definition = function(x, prob = 0.95){
            
            Mu = apply(x@mu,2,median)
            Delta = apply(x@delta,2,median)
            Phi = apply(x@phi,2,median)
            S = apply(x@s,2,median)
            Nu = apply(x@nu,2,median)
            Theta = apply(x@theta,2,median)
            
            HPDMu = coda::HPDinterval(coda::mcmc(x@mu), prob=prob)
            HPDDelta = coda::HPDinterval(coda::mcmc(x@delta), prob=prob)
            HPDPhi = coda::HPDinterval(coda::mcmc(x@phi), prob=prob)
            HPDS = coda::HPDinterval(coda::mcmc(x@s), prob=prob)
            HPDNu = coda::HPDinterval(coda::mcmc(x@nu), prob=prob)
            HPDTheta = coda::HPDinterval(coda::mcmc(x@theta), prob=prob)
            
            Output <- new("BASiCS_Summary", 
                          mu = cbind(Mu, HPDMu),
                          delta = cbind(Delta, HPDDelta),
                          phi = cbind(Phi, HPDPhi),
                          s = cbind(S, HPDS),
                          nu = cbind(Nu, HPDNu),
                          theta = cbind(Theta, HPDTheta))
            return(Output)
          })

#' @name plot-BASiCS_Chain-method
#' @aliases plot plot,BASiCS_Chain-method plot,BASiCS_Chain,ANY-method
#' 
#' @docType methods
#' @rdname plot-BASiCS_Chain-method
#' 
#' @title 'plot' method for BASiCS_Chain objects
#' 
#' @description 'plot' method for BASiCS_Chain objects
#' 
#' @param x A \code{BASiCS_Chain} object.
#' @param Param Name of the slot to be used for the plot. Possible values: \code{"mu"}, \code{"delta"}, 
#' \code{"phi"}, \code{"s"}, \code{"nu"} and \code{"theta"}
#' @param Gene Specifies which gene is requested. Required only if \code{Param = "mu"} or \code{"delta"} 
#' @param Cell Specifies which cell is requested. Required only if \code{Param = "phi", "s"} or \code{"nu"}
#' @param Batch Specifies which batch is requested. Required only if \code{Param = "theta"}
#' @param ylab As in \code{\link[graphics]{par}}. 
#' @param xlab As in \code{\link[graphics]{par}}. 
#' @param ... Other graphical parameters (see \code{\link[graphics]{par}}).
#' 
#' @examples
#' 
#' # See
#' help(BASiCS_MCMC)
#' 
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' 
#' @references 
#' 
#' Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. PLoS Computational Biology.
setMethod("plot",
          signature = "BASiCS_Chain",
          definition = function(x, 
                                Param = "mu",
                                Gene = NULL,
                                Cell = NULL, 
                                Batch = 1,
                                ylab = "",
                                xlab = "",
                                ...
          ){
            
            if(!(Param %in% c("mu", "delta", "phi", "s", "nu", "theta"))) stop("'Param' argument is invalid")
            if(Param %in% c("mu", "delta") & is.null(Gene))  stop("'Gene' value is required")
            if(Param %in% c("phi", "s", "nu") & is.null(Cell))  stop("'Cell' value is required")
            if(Param %in% c("theta") & is.null(Batch))  stop("'Batch' value is required")
            
            xlab = ifelse(xlab == "", "Iteration", xlab)
            
            if(Param == "mu") {object = x@mu; Column = Gene; if(ylab == "") ylab = bquote(mu[.(Column)])}
            if(Param == "delta") {object = x@delta; Column = Gene; if(ylab == "") ylab = bquote(delta[.(Column)])}
            if(Param == "phi") {object = x@phi; Column = Cell; if(ylab == "") ylab = bquote(phi[.(Column)])}
            if(Param == "s") {object = x@s; Column = Cell; ylab = if(ylab == "") ylab = bquote(s[.(Column)])}
            if(Param == "nu") {object = x@nu; Column = Cell; ylab = if(ylab == "") ylab = bquote(nu[.(Column)])} 
            if(Param == "theta") {object = x@theta; Column = Batch; ylab = if(ylab == "") ylab = bquote(theta[.(Column)])} 
            
            par(mfrow = c(1,2))
            plot(object[,Column], type="l", xlab = xlab, ylab = ylab, 
                 main = colnames(object)[Column], ...)
            stats::acf(object[,Column], main = "Autocorrelation")
          })


#' @name displayChainBASiCS-BASiCS_Chain-method
#' @aliases displayChainBASiCS displayChainBASiCS,BASiCS_Chain-method
#' 
#' @docType methods
#' @rdname displayChainBASiCS-BASiCS_Chain-method
#' 
#' @title Accessors for the slots of a BASiCS_Chain object
#' 
#' @description Accessors for the slots of a \code{\link[BASiCS]{BASiCS_Chain-class}}
#' 
#' @param object an object of class \code{\link[BASiCS]{BASiCS_Chain-class}}
#' @param Param Name of the slot to be used for the accessed. Possible values: \code{"mu"}, 
#' \code{"delta"}, \code{"phi"}, \code{"s"}, \code{"nu"} and \code{"theta"}.
#' 
#' @return The requested slot of an object of class \code{\link[BASiCS]{BASiCS_Chain-class}}
#' 
#' @examples
#' 
#' # See
#' help(BASiCS_MCMC)
#'   
#' @seealso \code{\link[BASiCS]{BASiCS_Chain-class}}
#' 
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' 
#' @references 
#' 
#' Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. PLoS Computational Biology.
setMethod("displayChainBASiCS",
          signature = "BASiCS_Chain",
          definition = function(object, 
                                Param = "mu"){
            
            if(!(Param %in% c("mu", "delta", "phi", "s", "nu", "theta"))) stop("'Param' argument is invalid")
            
            if(Param == "mu") {return(object@mu)}  
            if(Param == "delta") {return(object@delta)}  
            if(Param == "phi") {return(object@phi)}  
            if(Param == "s") {return(object@s)}  
            if(Param == "nu") {return(object@nu)}
            if(Param == "theta") {return(object@theta)} 
          })

##########################################################################
# Methods for BASiCS_Summary objects
##########################################################################

#' @name BASiCS_Summary-methods
#' @aliases show,BASiCS_Summary-method
#' 
#' @title 'show' method for BASiCS_Summary objects
#' 
#' @description 'show' method for \code{\link[BASiCS]{BASiCS_Summary-class}} objects.
#' 
#' @param object A \code{BASiCS_Summary} object.
#' 
#' @return Prints a summary of the properties of \code{object}.
#' 
#' @examples
#' 
#' # See
#' help(BASiCS_MCMC)
#' 
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' 
#' @references 
#' 
#' Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. PLoS Computational Biology.
#'
#' @rdname BASiCS_Summary-methods
setMethod("show",
          signature = "BASiCS_Summary",
          definition = function(object){
            q = nrow(object@mu)
            q.bio = nrow(object@delta)
            n = nrow(object@phi)
            nBatch = nrow(object@theta)
            message("An object of class ", class(object), "\n", sep = "")
            message(" Contains posterior medians and limits of HPD 95% interval for BASiCS parameters.\n")
            if(nBatch > 1)
            {
              message(" Dataset contains ", q, " genes (", q.bio, " biological and ", q-q.bio, " technical) and ", n, " cells (",nBatch," batches). \n", sep="")
            }
            else{message(" Dataset contains ", q, " genes (", q.bio, " biological and ", q-q.bio, " technical) and ", n, " cells (1 batch). \n", sep="")}
            message(" Elements (slots): mu, delta, phi, s, nu and theta.\n")
          })


#' @name plot-BASiCS_Summary-method
#' @aliases plot,BASiCS_Summary-method, plot,BASiCS_Summary,ANY-method
#' 
#' @docType methods
#' @rdname plot-BASiCS_Summary-method
#' 
#' @title 'plot' method for BASiCS_Summary objects
#' 
#' @description 'plot' method for BASiCS_Summary objects
#' 
#' @param x A \code{BASiCS_Summary} object.
#' @param Param Name of the slot to be used for the plot. Possible values: \code{"mu"}, \code{"delta"}, 
#' \code{"phi"}, \code{"s"}, \code{"nu"} and \code{"theta"}
#' @param Param2 Name of the second slot to be used for the plot. Possible values: \code{"mu"}, \code{"delta"}, 
#' \code{"phi"}, \code{"s"} and \code{"nu"} (combinations between gene-specific and cell-specific parameters are not admitted)
#' @param Genes Specifies which genes are requested. Required only if \code{Param = "mu"} or \code{"delta"} 
#' @param Cells Specifies which cells are requested. Required only if \code{Param = "phi", "s"} or \code{"nu"}
#' @param Batches Specifies which batches are requested. Required only if \code{Param = "theta"}
#' @param xlab As in \code{\link[graphics]{par}}. 
#' @param ylab As in \code{\link[graphics]{par}}.
#' @param xlim As in \code{\link[graphics]{par}}.  
#' @param ylim As in \code{\link[graphics]{par}}. 
#' @param pch As in \code{\link[graphics]{par}}. 
#' @param col As in \code{\link[graphics]{par}}. 
#' @param bty As in \code{\link[graphics]{par}}. 
#' @param SmoothPlot Logical parameter. If \code{TRUE}, transparency will be added to the color of the dots. 
#' @param ... Other graphical parameters (see \code{\link[graphics]{par}}).
#' 
#' @examples
#' 
#' # See
#' help(BASiCS_MCMC)
#' 
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' 
#' @references 
#' 
#' Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. PLoS Computational Biology.
#'
setMethod("plot",
          signature = "BASiCS_Summary",
          definition = function(x, 
                                Param = "mu",
                                Param2 = NULL,
                                Genes = NULL,
                                Cells = NULL, 
                                Batches = NULL,
                                xlab = "",
                                ylab = "",
                                xlim = "",
                                ylim = "",
                                pch = 16, 
                                col = "blue",
                                bty = "n",
                                SmoothPlot = TRUE,
                                ...
          ){
            
            if(!(Param %in% c("mu", "delta", "phi", "s", "nu", "theta"))) stop("'Param' argument is invalid")
            
            q = nrow(x@mu)
            q.bio = nrow(x@delta)
            n = nrow(x@phi)
            nBatch = nrow(x@theta)
            
            if(is.null(Genes)) {Genes = 1:q.bio}
            if(is.null(Cells)) {Cells = 1:n}
            if(is.null(Batches)) {Batches = 1:nBatch}
            
            if(is.null(Param2))
            {           
              if(Param == "mu") {object = x@mu; Columns = Genes; if(ylab == "") ylab = expression(mu[i]); if(xlab == "") xlab = "Gene"}
              if(Param == "delta") {object = x@delta; Columns = Genes; if(ylab == "") ylab = expression(delta[i]); if(xlab == "") xlab = "Gene"}
              if(Param == "phi") {object = x@phi; Columns = Cells; if(ylab == "") ylab = expression(phi[j]); if(xlab == "") xlab = "Cell"}
              if(Param == "s") {object = x@s; Columns = Cells; ylab = if(ylab == "") ylab = expression(s[j]); if(xlab == "") xlab = "Cell"}
              if(Param == "nu") {object = x@nu; Columns = Cells; ylab = if(ylab == "") ylab = expression(nu[j]); if(xlab == "") xlab = "Cell"} 
              if(Param == "theta") {object = x@theta; Columns = Batches; ylab = if(ylab == "") ylab = expression(theta[b]); if(xlab == "") xlab = "Batch"}
              
              if(ylim == "") {ylim = c(min(object[Columns,2]),max(object[Columns,3]))}
              
              plot(Columns, object[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)
              for(Column in Columns) 
              {
                BarLength = ifelse(length(Columns) <=10, 0.1, 2/length(Columns))
                segments(x0 = Column, y0 = object[Column,2], y1 = object[Column,3], col = col, ...)
                segments(x0 = Column - BarLength, y0 = object[Column,2], x1 = Column + BarLength, col = col, ...)
                segments(x0 = Column - BarLength, y0 = object[Column,3], x1 = Column + BarLength, col = col, ...)
              }
            }
            
            else{
              args = list(...)
              ValidCombination = FALSE
              if(SmoothPlot)
              {
                col = grDevices::rgb(grDevices::col2rgb(col)[1], grDevices::col2rgb(col)[2], 
                                     grDevices::col2rgb(col)[3],50,maxColorValue=255) 
              }
              
              if(Param == "mu" & Param2 == "delta")
              {
                ValidCombination = TRUE
                Columns = Genes; if(ylab == "") ylab = expression(delta[i]); if(xlab == "") xlab = expression(mu[i])
                if(xlim == "") {xlim = c(min(x@mu[Columns,1]),max(x@mu[Columns,1]))}
                if(ylim == "") {ylim = c(min(x@delta[Columns,1]),max(x@delta[Columns,1]))}

                plot(x@mu[Columns,1], x@delta[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...) 
              }
              if(Param == "delta" & Param2 == "mu")
              {
                ValidCombination = TRUE
                Columns = Genes; if(ylab == "") ylab = expression(mu[i]); if(xlab == "") xlab = expression(delta[i])
                if(ylim == "") {ylim = c(min(x@mu[Columns,1]),max(x@mu[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@delta[Columns,1]),max(x@delta[Columns,1]))}
                plot(x@delta[Columns,1], x@mu[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)
              }
              if(Param == "phi" & Param2 == "s")
              {
                ValidCombination = TRUE
                Columns = Cells; if(ylab == "") ylab = expression(s[i]); if(xlab == "") xlab = expression(phi[i])
                if(xlim == "") {xlim = c(min(x@phi[Columns,1]),max(x@phi[Columns,1]))}
                if(ylim == "") {ylim = c(min(x@s[Columns,1]),max(x@s[Columns,1]))}
                plot(x@phi[Columns,1], x@s[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)                
                }
              if(Param == "s" & Param2 == "phi")
              {
                ValidCombination = TRUE
                Columns = Cells; if(ylab == "") ylab = expression(phi[i]); if(xlab == "") xlab = expression(s[i])
                if(ylim == "") {ylim = c(min(x@phi[Columns,1]),max(x@phi[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@s[Columns,1]),max(x@s[Columns,1]))}
                plot(x@s[Columns,1], x@phi[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)                
                }
              if(Param == "phi" & Param2 == "nu")
              {
                ValidCombination = TRUE
                Columns = Cells; if(ylab == "") ylab = expression(nu[i]); if(xlab == "") xlab = expression(phi[i])
                if(xlim == "") {xlim = c(min(x@phi[Columns,1]),max(x@phi[Columns,1]))}
                if(ylim == "") {ylim = c(min(x@nu[Columns,1]),max(x@nu[Columns,1]))}
                plot(x@phi[Columns,1], x@nu[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)                
              }
              if(Param == "nu" & Param2 == "phi")
              {
                ValidCombination = TRUE
                Columns = Cells; if(ylab == "") ylab = expression(phi[i]); if(xlab == "") xlab = expression(nu[i])
                if(ylim == "") {ylim = c(min(x@phi[Columns,1]),max(x@phi[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@nu[Columns,1]),max(x@nu[Columns,1]))}
                plot(x@nu[Columns,1], x@phi[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)                
              }
              if(Param == "s" & Param2 == "nu")
              {
                ValidCombination = TRUE
                Columns = Cells; if(ylab == "") ylab = expression(nu[i]); if(xlab == "") xlab = expression(s[i])
                if(xlim == "") {xlim = c(min(x@s[Columns,1]),max(x@s[Columns,1]))}
                if(ylim == "") {ylim = c(min(x@nu[Columns,1]),max(x@nu[Columns,1]))}

                plot(x@s[Columns,1], x@nu[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)                
              }
              if(Param == "nu" & Param2 == "s")
              {
                ValidCombination = TRUE
                Columns = Cells; if(ylab == "") ylab = expression(s[i]); if(xlab == "") xlab = expression(nu[i])
                if(ylim == "") {ylim = c(min(x@s[Columns,1]),max(x@s[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@nu[Columns,1]),max(x@nu[Columns,1]))}

                plot(x@nu[Columns,1], x@s[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)                
              }
              
              if(ValidCombination == FALSE) {stop("Invalid combination for Param and Param2 \n - 
                                              Combinations between gene-specific and cell-specific parameters are not admitted")}
              }
            })


#' @name displaySummaryBASiCS-BASiCS_Summary-method
#' @aliases displaySummaryBASiCS displaySummaryBASiCS,BASiCS_Summary-method
#' 
#' @docType methods
#' @rdname displaySummaryBASiCS-BASiCS_Summary-method
#' 
#' @title Accessors for the slots of a BASiCS_Summary object
#' 
#' @description Accessors for the slots of a \code{\link[BASiCS]{BASiCS_Summary-class}}
#' 
#' @param object an object of class \code{\link[BASiCS]{BASiCS_Summary-class}}
#' @param Param Name of the slot to be used for the accessed. Possible values: \code{"mu"}, 
#' \code{"delta"}, \code{"phi"}, \code{"s"}, \code{"nu"} and \code{"theta"}
#'  
#' @return The requested slot of an object of class \code{\link[BASiCS]{BASiCS_Summary-class}}
#' 
#' @examples
#' 
#' # See
#' help(BASiCS_MCMC)
#'   
#' @seealso \code{\link[BASiCS]{BASiCS_Summary-class}}
#' 
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' 
#' @references 
#' 
#' Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. PLoS Computational Biology.
setMethod("displaySummaryBASiCS",
          signature = "BASiCS_Summary",
          definition = function(object, 
                                Param = "mu"){
            
            if(!(Param %in% c("mu", "delta", "phi", "s", "nu", "theta"))) stop("'Param' argument is invalid")
            
            if(Param == "mu") {return(object@mu)}  
            if(Param == "delta") {return(object@delta)}  
            if(Param == "phi") {return(object@phi)}  
            if(Param == "s") {return(object@s)}  
            if(Param == "nu") {return(object@nu)}   
            if(Param == "theta") {return(object@theta)}
          })