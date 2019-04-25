HiddenFindRBFLocations <- function(Data, k = 12) {
  # x <- log(rowMeans(SingleCellExperiment::counts(Data)))
  x <- log(Data)
  ml <- numeric(k - 2)
  range <- max(x) - min(x)
  ml[[1]] <- min(x)
  # arma::vec ml = arma::zeros(k - 2);
  # ml(0) = x.min();

  for (i in 2:(k-2)) {
    ml[[i]] = ml[[i-1]] + range / (k - 3);
  }
  ml
}
