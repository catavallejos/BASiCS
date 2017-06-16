##########################################################################
# New generic functions
##########################################################################
setGeneric("displayChainBASiCS", function(object, ...){})
setGeneric("displaySummaryBASiCS", function(object, ...){})
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
#' MCMC_Output <- BASiCS_MCMC(Data, N = 50, Thin = 2, Burn = 2)
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
#' 
#' @rdname BASiCS_Chain-methods
setMethod("show",
          signature = "BASiCS_Chain",
          definition = function(object){
            N = nrow(object@mu)
            q.bio = ncol(object@delta)
            n = ncol(object@phi)
            nBatch = ncol(object@theta)
            cat("An object of class ", class(object), "\n", sep = "")
            cat(" ", N," MCMC samples.\n", sep = "")
            if(nBatch > 1)
            {
              cat(" Dataset contains ", q.bio, " biological genes and ", n, " cells (",nBatch," batches). \n", sep="")
            }
            else{cat(" Dataset contains ", q.bio, " biological genes and ", n, " cells (1 batch). \n", sep="")}
            cat(" Elements (slots): mu, delta, phi, s, nu and theta.\n")
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
#' Data = makeExampleBASiCS_Data()
#' MCMC_Output <- BASiCS_MCMC(Data, N = 50, Thin = 2, Burn = 2)
#' MCMC_Summary <- Summary(MCMC_Output)
#' 
#' # See documentation of function BASiCS_MCMC
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. 
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
#' @aliases plot plot,BASiCS_Chain-method
#' 
#' @docType methods
#' @rdname plot-BASiCS_Chain-method
#' 
#' @title 'plot' method for BASiCS_Chain objects
#' 
#' @description 'plot' method for BASiCS_Chain objects
#' 
#' @param x A \code{BASiCS_Chain} object.
#' @param Param Name of the slot to be used for the plot. Possible values: \code{mu, delta, phi, s, nu, theta}
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
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. 
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
            
            plot(object[,Column], type="l", xlab = xlab, ylab = ylab, ...)
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
#' @param Param Name of the slot to be used for the accessed. Possible values: \code{mu, delta, phi, s, nu, theta}
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
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
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
#' Data = makeExampleBASiCS_Data()
#' MCMC_Output <- BASiCS_MCMC(Data, N = 50, Thin = 2, Burn = 2)
#' Summary(MCMC_Output)
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
#' 
#' @rdname BASiCS_Summary-methods
setMethod("show",
          signature = "BASiCS_Summary",
          definition = function(object){
            q = nrow(object@mu)
            q.bio = nrow(object@delta)
            n = nrow(object@phi)
            nBatch = nrow(object@theta)
            cat("An object of class ", class(object), "\n", sep = "")
            cat(" Contains posterior medians and limits of HPD 95% interval for BASiCS parameters.\n")
            if(nBatch > 1)
            {
              cat(" Dataset contains ", q, " genes (", q.bio, " biological and ", q-q.bio, " technical) and ", n, " cells (",nBatch," batches). \n", sep="")
            }
            else{cat(" Dataset contains ", q, " genes (", q.bio, " biological and ", q-q.bio, " technical) and ", n, " cells (1 batch). \n", sep="")}
            cat(" Elements (slots): mu, delta, phi, s, nu and theta.\n")
          })


#' @name plot-BASiCS_Summary-method
#' @aliases plot,BASiCS_Summary-method
#' 
#' @docType methods
#' @rdname plot-BASiCS_Summary-method
#' 
#' @title 'plot' method for BASiCS_Summary objects
#' 
#' @description 'plot' method for BASiCS_Summary objects
#' 
#' @param x A \code{BASiCS_Summary} object.
#' @param Param Name of the slot to be used for the plot. Possible values: \code{mu, delta, phi, s, nu, theta}
#' @param Param2 Name of the second slot to be used for the plot. Possible values: \code{mu, delta, phi, s, nu} (combinations between gene-specific and cell-specific parameters are not admitted)
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
#' @param ... Other graphical parameters (see \code{\link[graphics]{par}}).
#' 
#' @examples
#' 
#' # See
#' help(BASiCS_MCMC)
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. 
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
              ValidCombination = F
              if(Param == "mu" & Param2 == "delta")
              {
                ValidCombination = T
                Columns = Genes; if(ylab == "") ylab = expression(delta[i]); if(xlab == "") xlab = expression(mu[i])
                if(xlim == "") {xlim = c(min(x@mu[Columns,1]),max(x@mu[Columns,1]))}
                if(ylim == "") {ylim = c(min(x@delta[Columns,1]),max(x@delta[Columns,1]))}
                plot(x@mu[Columns,1], x@delta[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)                
              }
              if(Param == "delta" & Param2 == "mu")
              {
                ValidCombination = T
                Columns = Genes; if(ylab == "") ylab = expression(mu[i]); if(xlab == "") xlab = expression(delta[i])
                if(ylim == "") {ylim = c(min(x@mu[Columns,1]),max(x@mu[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@delta[Columns,1]),max(x@delta[Columns,1]))}
                plot(x@delta[Columns,1], x@mu[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)
              }
              if(Param == "phi" & Param2 == "s")
              {
                ValidCombination = T
                Columns = Cells; if(ylab == "") ylab = expression(s[i]); if(xlab == "") xlab = expression(phi[i])
                if(xlim == "") {xlim = c(min(x@phi[Columns,1]),max(x@phi[Columns,1]))}
                if(ylim == "") {ylim = c(min(x@s[Columns,1]),max(x@s[Columns,1]))}
                plot(x@phi[Columns,1], x@s[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)                
              }
              if(Param == "s" & Param2 == "phi")
              {
                ValidCombination = T
                Columns = Cells; if(ylab == "") ylab = expression(phi[i]); if(xlab == "") xlab = expression(s[i])
                if(ylim == "") {ylim = c(min(x@phi[Columns,1]),max(x@phi[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@s[Columns,1]),max(x@s[Columns,1]))}
                plot(x@s[Columns,1], x@phi[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)                
              }
              if(Param == "phi" & Param2 == "nu")
              {
                ValidCombination = T
                Columns = Cells; if(ylab == "") ylab = expression(nu[i]); if(xlab == "") xlab = expression(phi[i])
                if(xlim == "") {xlim = c(min(x@phi[Columns,1]),max(x@phi[Columns,1]))}
                if(ylim == "") {ylim = c(min(x@nu[Columns,1]),max(x@nu[Columns,1]))}
                plot(x@phi[Columns,1], x@nu[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)                
              }
              if(Param == "nu" & Param2 == "phi")
              {
                ValidCombination = T
                Columns = Cells; if(ylab == "") ylab = expression(phi[i]); if(xlab == "") xlab = expression(nu[i])
                if(ylim == "") {ylim = c(min(x@phi[Columns,1]),max(x@phi[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@nu[Columns,1]),max(x@nu[Columns,1]))}
                plot(x@nu[Columns,1], x@phi[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)                
              }
              if(Param == "s" & Param2 == "nu")
              {
                ValidCombination = T
                Columns = Cells; if(ylab == "") ylab = expression(nu[i]); if(xlab == "") xlab = expression(s[i])
                if(xlim == "") {xlim = c(min(x@s[Columns,1]),max(x@s[Columns,1]))}
                if(ylim == "") {ylim = c(min(x@nu[Columns,1]),max(x@nu[Columns,1]))}
                plot(x@s[Columns,1], x@nu[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)                
              }
              if(Param == "nu" & Param2 == "s")
              {
                ValidCombination = T
                Columns = Cells; if(ylab == "") ylab = expression(s[i]); if(xlab == "") xlab = expression(nu[i])
                if(ylim == "") {ylim = c(min(x@s[Columns,1]),max(x@s[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@nu[Columns,1]),max(x@nu[Columns,1]))}
                plot(x@nu[Columns,1], x@s[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)                
              }
              
              if(ValidCombination == F) {stop("Invalid combination for Param and Param2 \n - 
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
#' @param Param Name of the slot to be used for the accessed. Possible values: \code{mu, delta, phi, s, nu, theta}
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
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
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

##########################################################################
# Methods for BASiCS_D_Chain objects
##########################################################################

#' @name BASiCS_D_Chain-methods
#' @aliases show,BASiCS_D_Chain-method
#' 
#' @title 'show' method for BASiCS_D_Chain objects
#' 
#' @description 'show' method for \code{\link[BASiCS]{BASiCS_D_Chain-class}} objects.
#' 
#' @param object A \code{BASiCS_D_Chain} object.
#' 
#' @return Prints a summary of the properties of \code{object}.
#' 
#' @examples
#' 
#' Data = makeExampleBASiCS_D_Data()
#' #MCMC_Output <- BASiCS_D_MCMC(Data, N = 50, Thin = 2, Burn = 2)
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' @rdname BASiCS_D_Chain-methods
setMethod("show",
          signature = "BASiCS_D_Chain",
          definition = function(object){
            N = nrow(object@muTest)
            q.bio = ncol(object@muTest)
            n = ncol(object@phi)
            cat("An object of class ", class(object), "\n", sep = "")
            cat(" ", N," MCMC samples.\n", sep = "")
            cat(" Dataset contains ", q.bio, " biological genes and ", n, " cells (in total across both samples).\n", sep="")
            cat(" Offset = ", object@offset, ".\n", sep = "")
            cat(" Elements (slots): muTest, muRef, deltaTest, deltaRef, phi, s, nu, thetaTest, thetaRef and offset.\n")
            
          })

#' @name plot-BASiCS_D_Chain-method
#' @aliases plot plot,BASiCS_D_Chain-method
#' 
#' @docType methods
#' @rdname plot-BASiCS_D_Chain-method
#' 
#' @title 'plot' method for BASiCS_D_Chain objects
#' 
#' @description  'plot' method for BASiCS_D_Chain objects
#' 
#' @param x A \code{BASiCS_D_Chain} object.
#' @param Param Name of the slot to be used for the plot. Possible values: \code{muTest, muRef, deltaTest, deltaRef, phi, s, nu, thetaTest, thetaRef}
#' @param Gene Specifies which gene is requested. Required only if \code{Param = "muTest"}, \code{Param = "muRef"}, \code{"deltaTest"} or \code{"deltaRef"}
#' @param Cell Specifies which cell is requested. Required only if \code{Param = "phi", "s"} or \code{"nu"}
#' @param ylab As in \code{\link[graphics]{par}}. 
#' @param xlab As in \code{\link[graphics]{par}}. 
#' @param main As in \code{\link[graphics]{par}}. 
#' @param ... Other graphical parameters (see \code{\link[graphics]{par}}).
#' 
#' @examples
#' 
#' # See
#' help(BASiCS_D_MCMC)
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
setMethod("plot",
          signature = "BASiCS_D_Chain",
          definition = function(x, 
                                Param = "muTest",
                                Gene = NULL,
                                Cell = NULL, 
                                ylab = "",
                                xlab = "",
                                main = "",
                                ...
          ){
            
            if(!(Param %in% c("muTest", "muRef", 
                              "deltaTest", "deltaRef", 
                              "phi", "s", "nu", 
                              "thetaTest", "thetaRef"))) stop("'Param' argument is invalid")
            if(Param %in% c("muTest", "muRef", "deltaTest", "deltaRef") & is.null(Gene))  stop("'Gene' value is required")
            if(Param %in% c("phi", "s", "nu") & is.null(Cell))  stop("'Cell' value is required")
            
            xlab = ifelse(xlab == "", "Iteration", xlab)
            
            if(Param == "thetaTest") 
            {
              ylab = ifelse(ylab == "", expression(theta[test]), ylab)
              if(main == "") main = "Test group"
              plot(x@thetaTest, type="l", ylab = ylab, xlab = xlab, main = main, ...)
            }
            if(Param == "thetaRef") 
            {
              ylab = ifelse(ylab == "", expression(theta[ref]), ylab)
              if(main == "") main = "Reference group"
              plot(x@thetaRef, type="l", ylab = ylab, xlab = xlab, main = main, ...)
            }
            if(Param %in% c("muTest", "muRef", "deltaTest", "deltaRef", "phi", "s", "nu"))
            {
              if(Param == "muTest")    {object = x@muTest; Column = Gene; if(ylab == "") ylab = bquote(mu[.(Column)]); if(main == "") main = "Test group"}
              if(Param == "muRef")     {object = x@muRef; Column = Gene; if(ylab == "") ylab = bquote(tau[.(Column)]); if(main == "") main = "Reference group"}
              if(Param == "deltaTest") {object = x@deltaTest; Column = Gene; if(ylab == "") ylab = bquote(delta[.(Column)]); if(main == "") main = "Test group"}
              if(Param == "deltaRef")  {object = x@deltaRef; Column = Gene; if(ylab == "") ylab = bquote(omega[.(Column)]); if(main == "") main = "Reference group"}
              if(Param == "phi")       {object = x@phi; Column = Cell; if(ylab == "") ylab = bquote(phi[.(Column)])}
              if(Param == "s")         {object = x@s; Column = Cell; ylab = if(ylab == "") ylab = bquote(s[.(Column)])}
              if(Param == "nu")        {object = x@nu; Column = Cell; ylab = if(ylab == "") ylab = bquote(nu[.(Column)])} 
              
              plot(object[,Column], type="l", xlab = xlab, ylab = ylab, main = main, ...)
            }
          })

#' @name Summary
#' @aliases Summary Summary,BASiCS_D_Chain-method
#' 
#' @docType methods
#' @rdname Summary-BASiCS_D_Chain-method
#' 
#' @title 'Summary' method for BASiCS_D_Chain objects
#' 
#' @description For each of the BASiCS parameters, 
#' \code{Summary} returns the corresponding postior medians and limits of the high posterior
#' density interval with probabilty equal to \code{prob}.
#' 
#' @param x A \code{BASiCS_D_Chain} object.
#' @param prob \code{prob} argument for \code{\link[coda]{HPDinterval}} function. 
#' 
#' @return An object of class \code{\link[BASiCS]{BASiCS_D_Summary-class}}. 
#' 
#' @examples 
#' 
#' # See help(BASiCS_D_MCMC)
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
setMethod("Summary",
          signature = "BASiCS_D_Chain",
          definition = function(x, prob = 0.95){
            
            MuTest = apply(x@muTest,2,median)
            MuRef = apply(x@muRef,2,median)
            DeltaTest = apply(x@deltaTest,2,median)
            DeltaRef = apply(x@deltaRef,2,median)
            Phi = apply(x@phi,2,median)
            S = apply(x@s,2,median)
            Nu = apply(x@nu,2,median)
            ThetaTest = apply(x@thetaTest,2,median)
            ThetaRef = apply(x@thetaRef,2,median)
            
            HPDMuTest = coda::HPDinterval(coda::mcmc(x@muTest), prob=prob)
            HPDMuRef = coda::HPDinterval(coda::mcmc(x@muRef), prob=prob)
            HPDDeltaTest = coda::HPDinterval(coda::mcmc(x@deltaTest), prob=prob)
            HPDDeltaRef = coda::HPDinterval(coda::mcmc(x@deltaRef), prob=prob)
            HPDPhi = coda::HPDinterval(coda::mcmc(x@phi), prob=prob)
            HPDS = coda::HPDinterval(coda::mcmc(x@s), prob=prob)
            HPDNu = coda::HPDinterval(coda::mcmc(x@nu), prob=prob)
            HPDThetaTest = coda::HPDinterval(coda::mcmc(x@thetaTest), prob=prob)
            HPDThetaRef = coda::HPDinterval(coda::mcmc(x@thetaRef), prob=prob)
            
            Output <- new("BASiCS_D_Summary", 
                          muTest = cbind(MuTest, HPDMuTest),
                          muRef = cbind(MuRef, HPDMuRef),
                          deltaTest = cbind(DeltaTest, HPDDeltaTest),
                          deltaRef = cbind(DeltaRef, HPDDeltaRef),
                          phi = cbind(Phi, HPDPhi),
                          s = cbind(S, HPDS),
                          nu = cbind(Nu, HPDNu),
                          thetaTest = cbind(ThetaTest, HPDThetaTest),
                          thetaRef = cbind(ThetaRef, HPDThetaRef),
                          offset = x@offset,
                          probHPD = prob)
            show(Output)
            return(Output)
          })


#' @name displayChainBASiCS-BASiCS_D_Chain-method
#' @aliases displayChainBASiCS displayChainBASiCS,BASiCS_D_Chain-method
#' 
#' @docType methods
#' @rdname displayChainBASiCS-BASiCS_D_Chain-method
#' 
#' @title Accessors for the slots of a BASiCS_D_Chain object
#' 
#' @description Accessors for the slots of a \code{\link[BASiCS]{BASiCS_D_Chain-class}}
#' 
#' @param object an object of class \code{\link[BASiCS]{BASiCS_D_Chain-class}}
#' @param Param Name of the slot to be used for the accessed. 
#' Possible values: \code{muTest, muRef, deltaTest, deltaRef, phi, s, nu, thetaTest, thetaRef}
#' 
#' @return The requested slot of an object of class \code{\link[BASiCS]{BASiCS_D_Chain-class}}
#' 
#' @examples
#' 
#' # See
#' help(BASiCS_D_MCMC)
#'   
#' @seealso \code{\link[BASiCS]{BASiCS_D_Chain-class}}
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
setMethod("displayChainBASiCS",
          signature = "BASiCS_D_Chain",
          definition = function(object, 
                                Param = "muTest"){
            
            if(!(Param %in% c("muTest", "muRef", "deltaTest", "deltaRef", 
                              "phi", "s", "nu", 
                              "thetaTest", "thetaRef"))) stop("'Param' argument is invalid")
            
            if(Param == "muTest") {return(object@muTest)}  
            if(Param == "muRef") {return(object@muRef)} 
            if(Param == "deltaTest") {return(object@deltaTest)} 
            if(Param == "deltaRef") {return(object@deltaRef)}
            if(Param == "phi") {return(object@phi)}  
            if(Param == "s") {return(object@s)}  
            if(Param == "nu") {return(object@nu)}
            if(Param == "thetaTest") {return(object@thetaTest)} 
            if(Param == "thetaRef") {return(object@thetaRef)} 
          })

##########################################################################
# Methods for BASiCS_D_Summary objects
##########################################################################

#' @name BASiCS_D_Summary-methods
#' @aliases show,BASiCS_D_Summary-method
#' 
#' @title 'show' method for BASiCS_D_Summary objects
#' 
#' @description 'show' method for \code{\link[BASiCS]{BASiCS_D_Summary-class}} objects.
#' 
#' @param object A \code{BASiCS_D_Summary} object.
#' 
#' @return Prints a summary of the properties of \code{object}.
#' 
#' @examples
#' 
#' # see help(BASiCS_D_MCMC)
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @rdname BASiCS_D_Summary-methods
setMethod("show",
          signature = "BASiCS_D_Summary",
          definition = function(object){
            q.bio = nrow(object@muTest)
            q = nrow(object@deltaTest)
            n = nrow(object@phi)
            cat("An object of class ", class(object), "\n", sep = "")
            cat(" Contains posterior medians and limits of HPD 95% interval for BASiCSD parameters.\n")
            cat(" Dataset contains ", q, " genes (", q.bio, " biological and ", q-q.bio, " technical) and ", n, " cells.\n", sep="")
            cat(" Offset = ", object@offset, ".\n", sep = "")
            cat(" HPD interval probability = ", object@probHPD, ".\n", sep = "")
            cat(" Elements (slots): muTest, muRef, deltaTest, omegaRef, phi, s, nu, thetaTest, thetaRef, offset and probHPD.\n")
          })

#' @name plot-BASiCS_D_Summary-method
#' @aliases plot,BASiCS_D_Summary-method
#' 
#' @docType methods
#' @rdname plot-BASiCS_D_Summary-method
#' 
#' @title 'plot' method for BASiCS_D_Summary objects
#' 
#' @description 'plot' method for BASiCS_D_Summary objects
#' 
#' @param x A \code{BASiCS_D_Summary} object.
#' @param Param Name of the slot to be used for the plot. 
#' Possible values: \code{"muTest", "muRef", "deltaTest", "deltaRef", "phi", "s", "nu", "thetaTest", "thetaRef"}
#' @param Param2 Name of the second slot to be used for the plot. 
#' Possible values: \code{"muTest", "muRef", "deltaTest", "deltaRef", "phi", "s", "nu"} 
#' (combinations between gene-specific and cell-specific parameters are not allowed)
#' @param Genes Specifies which genes are requested. 
#' Required only if \code{Param = "muTest", "muRef", "deltaTest"} or \code{"deltaRef"} 
#' @param Cells Specifies which cells are requested. 
#' Required only if \code{Param = "phi", "s"} or \code{"nu"}
#' @param xlab As in \code{\link[graphics]{par}}. 
#' @param ylab As in \code{\link[graphics]{par}}.
#' @param xlim As in \code{\link[graphics]{par}}.  
#' @param ylim As in \code{\link[graphics]{par}}. 
#' @param main As in \code{\link[graphics]{par}}. 
#' @param pch As in \code{\link[graphics]{par}}. 
#' @param col As in \code{\link[graphics]{par}}. 
#' @param bty As in \code{\link[graphics]{par}}. 
#' @param ... Other graphical parameters (see \code{\link[graphics]{par}}).
#' 
#' @examples
#' 
#' # See
#' help(BASiCS_D_MCMC)
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
setMethod("plot",
          signature = "BASiCS_D_Summary",
          definition = function(x, 
                                Param = "muTest",
                                Param2 = NULL,
                                Genes = NULL,
                                Cells = NULL, 
                                xlab = "",
                                ylab = "",
                                xlim = "",
                                ylim = "",
                                main = "",
                                pch = 16, 
                                col = "blue",
                                bty = "n",
                                ...
          ){
            
            if(!(Param %in% c("muTest", "muRef", 
                              "deltaTest", "deltaRef", 
                              "phi", "s", "nu", 
                              "thetaTest", "thetaRef"))) stop("'Param' argument is invalid")
            
            q.bio = nrow(x@muTest)
            n = nrow(x@phi)
            
            if(is.null(Genes)) {Genes = 1:q.bio}
            if(is.null(Cells)) {Cells = 1:n}
            
            if(is.null(Param2))
            {           
              if(Param %in% c("thetaTest","thetaRef") ) 
              {  
                if(Param == "thetaTest") 
                {
                  ylab = ifelse(ylab == "", expression(theta[test]), ylab)
                  main = ifelse(main == "", "Test group", main)
                  if(ylim == "") {ylim = c(x@thetaTest[2],x@thetaTest[3])}
                  
                  plot(1, x@thetaTest[1], xlab = xlab, ylab = ylab, xlim = c(0,2), 
                       ylim = ylim, pch = pch, col = col, bty = bty, main = main, ...)
                  segments(x0 = 1, y0 = x@thetaTest[2], y1 = x@thetaTest[3], col = col, ...)
                  segments(x0 = 0.9, y0 = x@thetaTest[2], x1 = 1.1, col = col)
                  segments(x0 = 0.9, y0 = x@thetaTest[3], x1 = 1.1, col = col)                                    
                }
                if(Param == "thetaRef") 
                {
                  ylab = ifelse(ylab == "", expression(theta[ref]), ylab)
                  main = ifelse(main == "", "Reference group", main)
                  if(ylim == "") {ylim = c(x@thetaRef[2],x@thetaRef[3])}
                  
                  plot(1, x@thetaRef[1], xlab = xlab, ylab = ylab, xlim = c(0,2), 
                       ylim = ylim, pch = pch, col = col, bty = bty, main = main, ...)
                  segments(x0 = 1, y0 = x@thetaRef[2], y1 = x@thetaRef[3], col = col)
                  segments(x0 = 0.9, y0 = x@thetaRef[2], x1 = 1.1, col = col)
                  segments(x0 = 0.9, y0 = x@thetaRef[3], x1 = 1.1, col = col)                                    
                }
              }           
              else
              {
                if(Param == "muTest") {object = x@muTest; Columns = Genes; 
                                       if(ylab == "") ylab = expression(mu[i]); 
                                       if(xlab == "") xlab = "Gene"; 
                                       main = ifelse(main == "", "Test group", main)}
                if(Param == "muRef") {object = x@muRef; Columns = Genes; 
                                      if(ylab == "") ylab = expression(mu[i]); 
                                      if(xlab == "") xlab = "Gene"; 
                                      main = ifelse(main == "", "Reference group", main)}
                if(Param == "deltaTest") {object = x@deltaTest; Columns = Genes; 
                                          if(ylab == "") ylab = expression(delta[i]); 
                                          if(xlab == "") xlab = "Gene";
                                          main = ifelse(main == "", "Test group", main)}
                if(Param == "deltaRef") {object = x@deltaRef; Columns = Genes; 
                                         if(ylab == "") ylab = expression(delta[i]); 
                                         if(xlab == "") xlab = "Gene"; 
                                         main = ifelse(main == "", "Reference group", main)}
                if(Param == "phi") {object = x@phi; Columns = Cells; 
                                    if(ylab == "") ylab = expression(phi[j]); 
                                    if(xlab == "") xlab = "Cell"}
                if(Param == "s") {object = x@s; Columns = Cells; 
                                  ylab = if(ylab == "") ylab = expression(s[j]); 
                                  if(xlab == "") xlab = "Cell"}
                if(Param == "nu") {object = x@nu; Columns = Cells; 
                                   ylab = if(ylab == "") ylab = expression(nu[j]); 
                                   if(xlab == "") xlab = "Cell"} 
                
                if(ylim == "") {ylim = c(min(object[Columns,2]),max(object[Columns,3]))}
                
                plot(Columns, object[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main,...)
                for(Column in Columns) 
                {
                  segments(x0 = Column, y0 = object[Column,2], y1 = object[Column,3], col = col)
                  segments(x0 = Column-2/length(Columns), y0 = object[Column,2], x1 = Column+2/length(Columns), col = col)
                  segments(x0 = Column-2/length(Columns), y0 = object[Column,3], x1 = Column+2/length(Columns), col = col)
                }
              }
            }
            
            else{
              ValidCombination = F
              if(Param == "muTest" & Param2 == "deltaTest")
              {
                ValidCombination = T
                Columns = Genes; if(ylab == "") ylab = expression(delta[i]); if(xlab == "") xlab = expression(mu[i])
                if(ylim == "") {ylim = c(min(x@deltaTest[Columns,1]),max(x@deltaTest[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@muTest[Columns,1]),max(x@muTest[Columns,1]))}
                main = ifelse(main == "", "Test group", main)
                plot(x@muTest[Columns,1], x@deltaTest[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)
              }
              if(Param == "deltaTest" & Param2 == "muTest")
              {
                ValidCombination = T
                Columns = Genes; if(ylab == "") ylab = expression(mu[i]); if(xlab == "") xlab = expression(delta[i])
                if(ylim == "") {ylim = c(min(x@muTest[Columns,1]),max(x@muTest[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@deltaTest[Columns,1]),max(x@deltaTest[Columns,1]))}
                main = ifelse(main == "", "Test group", main)
                plot(x@deltaTest[Columns,1], x@muTest[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main,...)
              }
              
              if(Param == "muTest" & Param2 == "muRef")
              {
                ValidCombination = T
                Columns = Genes; if(ylab == "") ylab = expression(mu[i]); if(xlab == "") xlab = expression(mu[i])
                if(xlim == "") {xlim = c(min(x@muTest[Columns,1]),max(x@muTest[Columns,1]))}
                if(ylim == "") {ylim = c(min(x@muRef[Columns,1]),max(x@muRef[Columns,1]))}
                main = ifelse(main == "", "Test vs reference group", main)
                plot(x@muTest[Columns,1], x@muRef[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)                
              }
              if(Param == "muRef" & Param2 == "muTest")
              {
                ValidCombination = T
                Columns = Genes; if(ylab == "") ylab = expression(mu[i]); if(xlab == "") xlab = expression(mu[i])
                if(ylim == "") {ylim = c(min(x@muTest[Columns,1]),max(x@muTest[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@muRef[Columns,1]),max(x@muRef[Columns,1]))}
                main = ifelse(main == "", "Reference vs test group", main)
                plot(x@muRef[Columns,1], x@muTest[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)
              }
              
              if(Param == "deltaTest" & Param2 == "muRef")
              {
                ValidCombination = T
                Columns = Genes; if(ylab == "") ylab = expression(mu[i]); if(xlab == "") xlab = expression(delta[i])
                if(xlim == "") {xlim = c(min(x@deltaTest[Columns,1]),max(x@deltaTest[Columns,1]))}
                if(ylim == "") {ylim = c(min(x@muRef[Columns,1]),max(x@muRef[Columns,1]))}
                main = ifelse(main == "", "Test vs reference group", main)
                plot(x@deltaTest[Columns,1], x@muRef[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)                
              }
              if(Param == "muRef" & Param2 == "deltaTest")
              {
                ValidCombination = T
                Columns = Genes; if(ylab == "") ylab = expression(delta[i]); if(xlab == "") xlab = expression(mu[i])
                if(ylim == "") {ylim = c(min(x@deltaTest[Columns,1]),max(x@deltaTest[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@muRef[Columns,1]),max(x@muRef[Columns,1]))}
                main = ifelse(main == "", "Reference vs test group", main)
                plot(x@muRef[Columns,1], x@deltaTest[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)
              }
              
              if(Param == "muTest" & Param2 == "deltaRef")
              {
                ValidCombination = T
                Columns = Genes; if(ylab == "") ylab = expression(delta[i]); if(xlab == "") xlab = expression(mu[i])
                if(xlim == "") {xlim = c(min(x@muTest[Columns,1]),max(x@muTest[Columns,1]))}
                if(ylim == "") {ylim = c(min(x@deltaRef[Columns,1]),max(x@deltaRef[Columns,1]))}
                main = ifelse(main == "", "Test vs reference group", main)
                plot(x@muTest[Columns,1], x@deltaRef[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)                
              }
              if(Param == "deltaRef" & Param2 == "muTest")
              {
                ValidCombination = T
                Columns = Genes; if(ylab == "") ylab = expression(mu[i]); if(xlab == "") xlab = expression(delta[i])
                if(ylim == "") {ylim = c(min(x@muTest[Columns,1]),max(x@muTest[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@deltaRef[Columns,1]),max(x@deltaRef[Columns,1]))}
                main = ifelse(main == "", "Reference vs test group", main)
                plot(x@deltaRef[Columns,1], x@muTest[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)
              }
              
              if(Param == "deltaRef" & Param2 == "muRef")
              {
                ValidCombination = T
                Columns = Genes; if(ylab == "") ylab = expression(mu[i]); if(xlab == "") xlab = expression(delta[i])
                if(xlim == "") {xlim = c(min(x@deltaRef[Columns,1]),max(x@deltaRef[Columns,1]))}
                if(ylim == "") {ylim = c(min(x@muRef[Columns,1]),max(x@muRef[Columns,1]))}
                main = ifelse(main == "", "Reference group", main)
                plot(x@deltaRef[Columns,1], x@muRef[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)                
              }
              if(Param == "muRef" & Param2 == "deltaRef")
              {
                ValidCombination = T
                Columns = Genes; if(ylab == "") ylab = expression(delta[i]); if(xlab == "") xlab = expression(mu[i])
                if(ylim == "") {ylim = c(min(x@deltaRef[Columns,1]),max(x@deltaRef[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@muRef[Columns,1]),max(x@muRef[Columns,1]))}
                main = ifelse(main == "", "Reference group", main)
                plot(x@muRef[Columns,1], x@deltaRef[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)
              }
              
              if(Param == "deltaRef" & Param2 == "deltaTest")
              {
                ValidCombination = T
                Columns = Genes; if(ylab == "") ylab = expression(delta[i]); if(xlab == "") xlab = expression(delta[i])
                if(xlim == "") {xlim = c(min(x@deltaRef[Columns,1]),max(x@deltaRef[Columns,1]))}
                if(ylim == "") {ylim = c(min(x@deltaTest[Columns,1]),max(x@deltaTest[Columns,1]))}
                main = ifelse(main == "", "Reference vs test group", main)
                plot(x@deltaRef[Columns,1], x@deltaTest[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)                
              }
              if(Param == "deltaTest" & Param2 == "deltaRef")
              {
                ValidCombination = T
                Columns = Genes; if(ylab == "") ylab = expression(delta[i]); if(xlab == "") xlab = expression(delta[i])
                if(ylim == "") {ylim = c(min(x@deltaRef[Columns,1]),max(x@deltaRef[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@deltaTest[Columns,1]),max(x@deltaTest[Columns,1]))}
                main = ifelse(main == "", "Test vs reference group", main)
                plot(x@deltaTest[Columns,1], x@deltaRef[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)
              }
              
              if(Param == "phi" & Param2 == "s")
              {
                ValidCombination = T
                Columns = Cells; if(ylab == "") ylab = expression(s[i]); if(xlab == "") xlab = expression(phi[i])
                if(xlim == "") {xlim = c(min(x@phi[Columns,1]),max(x@phi[Columns,1]))}
                if(ylim == "") {ylim = c(min(x@s[Columns,1]),max(x@s[Columns,1]))}
                plot(x@phi[Columns,1], x@s[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)                
              }
              if(Param == "s" & Param2 == "phi")
              {
                ValidCombination = T
                Columns = Cells; if(ylab == "") ylab = expression(phi[i]); if(xlab == "") xlab = expression(s[i])
                if(ylim == "") {ylim = c(min(x@phi[Columns,1]),max(x@phi[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@s[Columns,1]),max(x@s[Columns,1]))}
                plot(x@s[Columns,1], x@phi[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)                
              }
              if(Param == "phi" & Param2 == "nu")
              {
                ValidCombination = T
                Columns = Cells; if(ylab == "") ylab = expression(nu[i]); if(xlab == "") xlab = expression(phi[i])
                if(xlim == "") {xlim = c(min(x@phi[Columns,1]),max(x@phi[Columns,1]))}
                if(ylim == "") {ylim = c(min(x@nu[Columns,1]),max(x@nu[Columns,1]))}
                plot(x@phi[Columns,1], x@nu[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)                
              }
              if(Param == "nu" & Param2 == "phi")
              {
                ValidCombination = T
                Columns = Cells; if(ylab == "") ylab = expression(phi[i]); if(xlab == "") xlab = expression(nu[i])
                if(ylim == "") {ylim = c(min(x@phi[Columns,1]),max(x@phi[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@nu[Columns,1]),max(x@nu[Columns,1]))}
                plot(x@nu[Columns,1], x@phi[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)                
              }
              if(Param == "s" & Param2 == "nu")
              {
                ValidCombination = T
                Columns = Cells; if(ylab == "") ylab = expression(nu[i]); if(xlab == "") xlab = expression(s[i])
                if(xlim == "") {xlim = c(min(x@s[Columns,1]),max(x@s[Columns,1]))}
                if(ylim == "") {ylim = c(min(x@nu[Columns,1]),max(x@nu[Columns,1]))}
                plot(x@s[Columns,1], x@nu[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)                
              }
              if(Param == "nu" & Param2 == "s")
              {
                ValidCombination = T
                Columns = Cells; if(ylab == "") ylab = expression(s[i]); if(xlab == "") xlab = expression(nu[i])
                if(ylim == "") {ylim = c(min(x@s[Columns,1]),max(x@s[Columns,1]))}
                if(xlim == "") {xlim = c(min(x@nu[Columns,1]),max(x@nu[Columns,1]))}
                plot(x@nu[Columns,1], x@s[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, 
                     pch = pch, col = col, bty = bty, main = main, ...)                
              }
              
              if(ValidCombination == F) {stop("Invalid combination for Param and Param2 \n - 
                                              Combinations between gene-specific and cell-specific parameters are not admitted")}
              }
            })

#' @name displaySummaryBASiCS-BASiCS_D_Summary-method
#' @aliases displaySummaryBASiCS displaySummaryBASiCS,BASiCS_D_Summary-method
#' 
#' @docType methods
#' @rdname displaySummaryBASiCS-BASiCS_D_Summary-method
#' 
#' @title Accessors for the slots of a BASiCS_D_Summary object
#' 
#' @description Accessors for the slots of a \code{\link[BASiCS]{BASiCS_D_Summary-class}}
#' 
#' @param object an object of class \code{\link[BASiCS]{BASiCS_D_Summary-class}}
#' @param Param Name of the slot to be used for the accessed. 
#' Possible values: \code{"muTest", "muRef", "deltaTest", "deltaRef", "phi", "s", "nu", "thetaTest", "thetaRef"}
#'  
#' @return The requested slot of an object of class \code{\link[BASiCS]{BASiCS_D_Summary-class}}
#' 
#' @examples
#' 
#' # See
#' help(BASiCS_D_MCMC)
#'   
#' @seealso \code{\link[BASiCS]{BASiCS_D_Summary-class}}
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
setMethod("displaySummaryBASiCS",
          signature = "BASiCS_D_Summary",
          definition = function(object, 
                                Param = "muTest"){
            
            if(!(Param %in% c("muTest", "muRef", 
                              "deltaTest", "deltaRef", 
                              "phi", "s", "nu", 
                              "thetaTest", "thetaRef"))) stop("'Param' argument is invalid")
            
            if(Param == "muTest") {return(object@muTest)}  
            if(Param == "muRef") {return(object@muRef)}  
            if(Param == "deltaTest") {return(object@deltaTest)}  
            if(Param == "deltaRef") {return(object@deltaRef)}
            if(Param == "phi") {return(object@phi)}  
            if(Param == "s") {return(object@s)}  
            if(Param == "nu") {return(object@nu)}
            if(Param == "thetaTest") {return(object@thetaTest)} 
            if(Param == "thetaRef") {return(object@thetaRef)} 
          })

