##########################################################################
# New generic functions
##########################################################################
setGeneric("displaySpikeInput", function(object){})
setGeneric("displayTechIndicator", function(object){})

##########################################################################
# Methods for BASiCS_Data objects
##########################################################################

#' @name BASiCS_Data-methods
#' @aliases show show,BASiCS_Data-method
#' 
#' @title S4 methods for BASiCS_Data objects
#' 
#' @description S4 methods for \code{\link[BASiCS]{BASiCS_Data-class}} objects.  
#' 
#' @param object A \code{BASiCS_Data} object.
#' @param type Only required for \code{counts} method. A string indicating which genes are returned. If \code{all}, returns expression counts for all genes. 
#' If \code{biological}, returns expression counts for biological genes only. If \code{technical}, returns expression 
#' counts for technical (spike-in) genes only. 
#' 
#' @return 
#' \describe{
#' \item{\code{show}}{Prints a summary of the properties of \code{object}.}
#' 
#' \item{\code{counts}}{If \code{type = "all"}, returns the \code{Counts} slot of \code{object}. 
#' If \code{type = "biological"} or \code{type = "spike-in"}, only biological or technical (spike-in) counts are returned, respectively.}
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
#' @seealso \code{\link[BASiCS]{BASiCS_Data-class}}, \code{\link[BASiCS]{BASiCS_MCMC_Start}}, \code{\link[BASiCS]{BASiCS_MCMC}}.
#' 
#' @rdname BASiCS_Data-methods
setMethod("show",
          signature = "BASiCS_Data",
          definition = function(object){
            q = nrow(object@Counts)
            n = ncol(object@Counts)
            q.bio = length(object@SpikeInput)
            cat("An object of class ", class(object), "\n", sep = "")
            cat(" Dataset contains ", q, " genes (", q.bio, " biological and ", q-q.bio, " technical) and ", n, " cells.\n", sep="")
            cat(" Elements (slots): Counts, Tech and SpikeInput.\n")
          })

#' @name BASiCS_Data-methods
#' @aliases counts counts,BASiCS_Data-method
#' 
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
#' @aliases Summary,BASiCS_Chain-method
#' 
#' @docType methods
#' @rdname Summary-methods
#' 
#' @title 'Summary' method for BASiCS_Chain objects
#' 
#' @description For each of the BASiCS parameters (see Vallejos et al 2015), 
#' \code{Summary} returns the corresponding postior medians and limits of the high posterior
#' density interval (probabilty equal to \code{prob}. )
#' 
#' @param x object A \code{BASiCS_Chain} object.
#' @param prob \code{prob} argument for function \code{\link[coda]{HPDinterval}}
#' 
#' @return An object of class \code{\link[BASiCS]{BASiCS_Summary-class}}. 
#' 
#' @examples 
#' 
#' # See documentation of function BASiCS_MCMC
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. 
setMethod("Summary",
          signature = "BASiCS_Chain",
          definition = function(x, prob = 0.95, ..., na.rm = FALSE){
            
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

#' @name plot
#' @aliases plot,BASiCS_Chain-method
#' 
#' @docType methods
#' @rdname plot-methods
#' 
#' @title 'plot' method for BASiCS_Chain and BASiCS_Summary objects
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. 
setMethod("plot",
          signature = "BASiCS_Chain",
          definition = function(x, y, 
                                Param = "mu",
                                Column = 1,
                                ylab = "",
                                xlab = "",
                                bty = "",
                                ...
                                ){
            
            if(!(Param %in% c("mu", "delta", "phi", "s", "nu", "theta"))) stop("'Param' argument is invalid")
            
            if(xlab == "") xlab = "Iteration"
            if(bty == "") bty = "n"
            
            if(Param == "theta") 
            {
                ylab = ifelse(ylab == "", expression(theta), ylab)
                plot(x@theta, type="l", bty = bty, ylab = ylab, xlab = xlab, ...)
            }
            else
            {
              if(Param == "mu") {object = x@mu; if(ylab == "") ylab = bquote(mu[.(Column)])}
              if(Param == "delta") {object = x@delta; if(ylab == "") ylab = bquote(delta[.(Column)])}
              if(Param == "phi") {object = x@phi; if(ylab == "") ylab = bquote(phi[.(Column)])}
              if(Param == "s") {object = x@s; ylab = if(ylab == "") ylab = bquote(s[.(Column)])}
              if(Param == "nu") {object = x@nu; ylab = if(ylab == "") ylab = bquote(nu[.(Column)])} 
                
              plot(object[,Column], type="l", bty=bty, xlab = xlab, ylab = ylab, ...)
            }
          })

##########################################################################
# Methods for BASiCS_Summary objects
##########################################################################

setMethod("show",
          signature = "BASiCS_Summary",
          definition = function(object){
            q = nrow(object@mu)
            q.bio = nrow(object@delta)
            n = nrow(object@phi)
            cat("An object of class ", class(object), "\n", sep = "")
            cat(" Contains posterior medians and limits of HPD 95% interval for BASiCS parameters")
            cat(" Dataset contains ", q, " genes (", q.bio, " biological and ", q-q.bio, " technical) and ", n, " cells.\n", sep="")
            cat(" Elements (slots): mu, delta, phi, s, nu and theta.\n")
          })

