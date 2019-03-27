#########################################################################
# New generic functions
setGeneric("displayChainBASiCS", 
           function(object, ...) standardGeneric("displayChainBASiCS"))
setGeneric("displaySummaryBASiCS", 
           function(object, ...) standardGeneric("displaySummaryBASiCS"))
setGeneric("BASiCS_showFit", function(object, ...) 
  standardGeneric("BASiCS_showFit"))

#if (!isGeneric("plot")) {
#  setGeneric("plot")
#  setMethod("plot", signature = "ANY", graphics::plot)
#}
