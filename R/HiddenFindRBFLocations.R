HiddenFindRBFLocations <- function(Data, k = 12) {
  as.vector(estimateRBFLocations(k, log(Data)))
}
