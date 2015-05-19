##########################################################################
# New generic functions
##########################################################################
setGeneric("displaySpikeInput", function(object){})
setGeneric("displayTechIndicator", function(object){})

##########################################################################
# Methods for BASiCS_Data objects
##########################################################################

#' @name BASiCS_Data-methods
#' @aliases show,BASiCS_Data-method
#' 
#' @title S4 methods for BASiCS_Data objects
#' 
#' @description S4 methods for \code{\link[BASiCS]{BASiCS_Data-class}} objects.  
#' 
#' @param object A \code{BASiCS_Data} object.
#' @param type Only required for \code{counts} method. A string indicating which genes are returned. 
#' 
#' @return 
#' \describe{
#' \item{\code{show}}{Prints a summary of the properties of \code{object}.}
#' 
#' \item{\code{counts}}{\describe{
#'  \item{\code{if(type = "all")}}{Returns the \code{Counts} slot of \code{object}.}
#'  \item{\code{if(type = "biological")}}{Returns the \code{Counts} slot of \code{object}, biological genes only.}
#'  \item{\code{if(type = "technical")}}{Returns the \code{Counts} slot of \code{object}, tehcnical genes only.}}}
#' 
#' \item{\code{displayTechIndicator}}{Returns \code{Tech} slot of \code{object}.}
#' 
#' \item{\code{displaySpikeInput}}{Returns \code{SpikeInput} slot of \code{object}.}
#' }
#' 
#' @examples
#' 
#' Data = makeExampleBASiCS_Data()
#' show(Data)
#' head(counts(Data))
#' dim(counts(Data, type="biological"))
#' dim(counts(Data, type="technical"))
#' displayTechIndicator(Data)
#' displaySpikeInput(Data)
#' 
#' @seealso \code{\link[BASiCS]{BASiCS_Data-class}}
#'  
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
#
#' @rdname BASiCS_Data-methods
setMethod("show",
          signature = "BASiCS_Data",
          definition = function(object){
            q = nrow(object@Counts)
            n = ncol(object@Counts)
            q.bio = nrow(object@Counts) - length(object@SpikeInput)
            cat("An object of class ", class(object), "\n", sep = "")
            cat(" Dataset contains ", q, " genes (", q.bio, " biological and ", q-q.bio, " technical) and ", n, " cells.\n", sep="")
            cat(" Elements (slots): Counts, Tech and SpikeInput.\n")
          })

#' @name BASiCS_Data-methods
#' @aliases counts counts,BASiCS_Data-method
#' @rdname BASiCS_Data-methods
setMethod("counts",
          signature = "BASiCS_Data",
          definition = function(object, type = "all"){
            
            if(!(type %in% c("all", "biological", "technical"))) stop("Invalid option for argument 'type'.")
            
            if(type == "all") return(object@Counts)
            if(type == "biological") return(object@Counts[!object@Tech,])
            if(type == "technical") return(object@Counts[object@Tech,])
          })

#' @name BASiCS_Data-methods
#' @aliases displaySpikeInput displaySpikeInput,BASiCS_Data-method
#' @rdname BASiCS_Data-methods
setMethod("displaySpikeInput",
          signature = "BASiCS_Data",
          definition = function(object){
                  return(object@SpikeInput)
          })

#' @name BASiCS_Data-methods
#' @aliases displayTechIndicator displayTechIndicator,BASiCS_Data-method
#' @rdname BASiCS_Data-methods
setMethod("displayTechIndicator",
          signature = "BASiCS_Data",
          definition = function(object){
            return(object@Tech)
          })

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
            q = ncol(object@mu)
            q.bio = ncol(object@delta)
            n = ncol(object@phi)
            cat("An object of class ", class(object), "\n", sep = "")
            cat(" ", N," MCMC samples.\n", sep = "")
            cat(" Dataset contains ", q, " genes (", q.bio, " biological and ", q-q.bio, " technical) and ", n, " cells.\n", sep="")
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
            Theta = median(x@theta)
            
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
#' @param x A \code{BASiCS_Chain} object.
#' @param Param Name of the slot to be used for the plot. Possible values: \code{mu, delta, phi, s, nu, theta}
#' @param Gene Specifies which gene is requested. Required only if \code{Param = "mu"} or \code{"delta"} 
#' @param Cell Specifies which cell is requested. Required only if \code{Param = "phi", "s"} or \code{"nu"}
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
                                ylab = "",
                                xlab = "",
                                ...
                                ){
            
            if(!(Param %in% c("mu", "delta", "phi", "s", "nu", "theta"))) stop("'Param' argument is invalid")
            if(Param %in% c("mu", "delta") & is.null(Gene))  stop("'Gene' value is required")
            if(Param %in% c("phi", "s", "nu") & is.null(Cell))  stop("'Cell' value is required")
            
            xlab = ifelse(xlab == "", "Iteration", xlab)
            
            if(Param == "theta") 
            {
                ylab = ifelse(ylab == "", expression(theta), ylab)
                plot(x@theta, type="l", ylab = ylab, xlab = xlab, ...)
            }
            else
            {
              if(Param == "mu") {object = x@mu; Column = Gene; if(ylab == "") ylab = bquote(mu[.(Column)])}
              if(Param == "delta") {object = x@delta; Column = Gene; if(ylab == "") ylab = bquote(delta[.(Column)])}
              if(Param == "phi") {object = x@phi; Column = Cell; if(ylab == "") ylab = bquote(phi[.(Column)])}
              if(Param == "s") {object = x@s; Column = Cell; ylab = if(ylab == "") ylab = bquote(s[.(Column)])}
              if(Param == "nu") {object = x@nu; Column = Cell; ylab = if(ylab == "") ylab = bquote(nu[.(Column)])} 
                
              plot(object[,Column], type="l", xlab = xlab, ylab = ylab, ...)
            }
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
            cat("An object of class ", class(object), "\n", sep = "")
            cat(" Contains posterior medians and limits of HPD 95% interval for BASiCS parameters.\n")
            cat(" Dataset contains ", q, " genes (", q.bio, " biological and ", q-q.bio, " technical) and ", n, " cells.\n", sep="")
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
#' @param x A \code{BASiCS_Summary} object.
#' @param Param Name of the slot to be used for the plot. Possible values: \code{mu, delta, phi, s, nu, theta}
#' @param Param2 Name of the second slot to be used for the plot. Possible values: \code{mu, delta, phi, s, nu} (combinations between gene-specific and cell-specific parameters are not admitted)
#' @param Genes Specifies which genes are requested. Required only if \code{Param = "mu"} or \code{"delta"} 
#' @param Cells Specifies which cells are requested. Required only if \code{Param = "phi", "s"} or \code{"nu"}
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
            
            if(is.null(Genes)) {Genes = 1:q.bio}
            if(is.null(Cells)) {Cells = 1:n}
            
            if(is.null(Param2))
            {           
              if(Param == "theta") 
              {              
                ylab = ifelse(ylab == "", expression(theta), ylab)
                if(ylim == "") {ylim = c(x@theta[2],x@theta[3])}
             
                plot(1, x@theta[1], xlab = xlab, ylab = ylab, xlim = c(0,2), ylim = ylim, pch = pch, col = col, bty = bty, ...)
                segments(x0 = 1, y0 = x@theta[2], y1 = x@theta[3], col = col, ...)
                segments(x0 = 0.9, y0 = x@theta[2], x1 = 1.1, col = col, ...)
                segments(x0 = 0.9, y0 = x@theta[3], x1 = 1.1, col = col, ...)
              }           
              else
              {
                if(Param == "mu") {object = x@mu; Columns = Genes; if(ylab == "") ylab = expression(mu[i]); if(xlab == "") xlab = "Gene"}
                if(Param == "delta") {object = x@delta; Columns = Genes; if(ylab == "") ylab = expression(delta[i]); if(xlab == "") xlab = "Gene"}
                if(Param == "phi") {object = x@phi; Columns = Cells; if(ylab == "") ylab = expression(phi[j]); if(xlab == "") xlab = "Cell"}
                if(Param == "s") {object = x@s; Columns = Cells; ylab = if(ylab == "") ylab = expression(s[j]); if(xlab == "") xlab = "Cell"}
                if(Param == "nu") {object = x@nu; Columns = Cells; ylab = if(ylab == "") ylab = expression(nu[j]); if(xlab == "") xlab = "Cell"} 
              
                if(ylim == "") {ylim = c(min(object[Columns,2]),max(object[Columns,3]))}
              
                plot(Columns, object[Columns,1], xlab = xlab, ylab = ylab, ylim = ylim, pch = pch, col = col, bty = bty, ...)
                for(Column in Columns) 
                {
                  segments(x0 = Column, y0 = object[Column,2], y1 = object[Column,3], col = col)
                  segments(x0 = Column-2/length(Columns), y0 = object[Column,2], x1 = Column+2/length(Columns), col = col, ...)
                  segments(x0 = Column-2/length(Columns), y0 = object[Column,3], x1 = Column+2/length(Columns), col = col, ...)
                }
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

