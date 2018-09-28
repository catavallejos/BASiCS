#' @importFrom stats acf
#' @importFrom Rcpp evalCpp
#' @importFrom coda HPDinterval mcmc
#' @importFrom methods Summary is new show slotNames .hasSlot
#' @importFrom testthat context test_check
#' @importFrom graphics abline barplot legend lines par points segments smoothScatter plot
#' @importFrom data.table fread
#' @importFrom S4Vectors metadata DataFrame
#' @importFrom SummarizedExperiment assay colData 
#' @importFrom SingleCellExperiment isSpike
#' @importFrom stats median model.matrix rgamma rpois runif acf
#' @importFrom utils packageVersion write.table
#' @importFrom matrixStats rowMeans2 colMeans2 rowVars colVars colMedians rowMedians
#' @importFrom grDevices colorRampPalette
#' @importFrom MASS mvrnorm
#' @importMethodsFrom BiocGenerics counts updateObject subset colnames rownames
#' @importMethodsFrom scran computeSumFactors
#' @useDynLib BASiCS
#' @exportPattern "^[^\\Hidden]"
#' @importFrom ggplot2 aes_string scale_fill_gradientn geom_point theme_minimal geom_line geom_ribbon aes aes_string geom_hex geom_line geom_point geom_ribbon ggplot
#' @import ggplot2
#' @import SingleCellExperiment
#' @import KernSmooth
NULL
