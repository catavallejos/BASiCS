#########################################################################
# New generic functions
setGeneric("displayChainBASiCS", 
           function(object, ...) standardGeneric("displayChainBASiCS"))
setGeneric("displaySummaryBASiCS", 
           function(object, ...) standardGeneric("displaySummaryBASiCS"))
setGeneric("BASiCS_showFit", function(object, ...) 
  standardGeneric("BASiCS_showFit"))
setGeneric("BASiCS_diagHist", function(object, ...) {
  standardGeneric("BASiCS_diagHist")
})
setGeneric("BASiCS_diagPlot", function(object, ...) {
  standardGeneric("BASiCS_diagPlot")
})

if (!isGeneric("plot")) {
#' @export
setGeneric("plot")
setMethod("plot", signature = "ANY", graphics::plot)
}
