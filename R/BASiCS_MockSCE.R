#' Create a mock \linkS4class{SingleCellExperiment} object.
#' 
#' Creates a \linkS4class{SingleCellExperiment} object of Poisson-distributed
#' approximating a homogeneous cell population.
#' 
#' @param NGenes Integer value specifying the number of genes that will be
#' present in the output.
#' @param NCells Integer value specifying the number of cells that will be
#' present in the output.
#' @param NSpikes Integer value specifying the number of spike-in genes that
#' will be present in the output.
#' @param WithBatch Logical value specifying whether a dummy \code{BatchInfo}
#' is included in the output.
#' @return A \linkS4class{SingleCellExperiment} object.
#' @examples
#' BASiCS_MockSCE()
#' @export
BASiCS_MockSCE <- function(NGenes = 100, NCells = 100, NSpikes = 20, WithBatch = TRUE) {
  mu <- rlnorm(NGenes)
  nu <- rgamma(NCells, 5, 5)
  mat <- mu %*% t(nu)
  counts <- matrix(
    rpois(length(mat), lambda = mat),
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = list(
      paste("Gene", seq_len(nrow(mat))),
      paste("Cell", seq_len(ncol(mat)))
    )
  )
  Data <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts)
  )
  if (WithBatch) {
    colData(Data)$BatchInfo <- sample(2, ncol(Data), replace = TRUE)
  }
  if (NSpikes > 0) {
    spikes <- rlnorm(NSpikes)
    spike_mat <- spikes %*% t(nu)
    spike_counts <- matrix(
      rpois(length(spike_mat), lambda = spike_mat),
      nrow = nrow(spike_mat),
      ncol = ncol(spike_mat),
      dimnames = list(
        paste("Spike", seq_len(nrow(spike_mat))),
        paste("Cell", seq_len(ncol(spike_mat)))
      )
    )
    Spikes <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = spike_counts)
    )
    rowData(Spikes)[[1]] <- rownames(Spikes)
    rowData(Spikes)[[2]] <- round(spikes) + 1
    altExp(Data, "spike-ins") <- Spikes
  }
  Data
}
