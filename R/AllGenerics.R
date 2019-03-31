#########################################################################
# New generic functions
setGeneric("displayChainBASiCS", 
           function(object, ...) standardGeneric("displayChainBASiCS"))
setGeneric("displaySummaryBASiCS", 
           function(object, ...) standardGeneric("displaySummaryBASiCS"))
if (!isGeneric("plot")) { setGeneric("plot") }