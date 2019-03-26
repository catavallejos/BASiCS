#' @useDynLib BASiCS
#' @exportPattern "^[^\\Hidden]"
#' @importFrom coda HPDinterval mcmc effectiveSize
#' @importFrom cowplot plot_grid
#' @importFrom data.table fread
#' @importFrom ggplot2 aes_string scale_fill_gradientn geom_point theme_minimal geom_line geom_ribbon aes aes_string geom_hex geom_line geom_point geom_ribbon ggplot
#' @importFrom graphics abline barplot legend lines par points segments smoothScatter plot
#' @importFrom grDevices colorRampPalette
#' @importFrom ggExtra ggMarginal
#' @importFrom MASS mvrnorm
#' @importFrom matrixStats rowMeans2 colMeans2 rowVars colVars colMedians rowMedians
#' @importFrom methods Summary is new show slotNames .hasSlot
#' @importFrom Rcpp evalCpp
#' @importFrom S4Vectors metadata DataFrame
#' @importFrom stats median model.matrix rgamma rpois runif acf
#' @importFrom SingleCellExperiment isSpike
#' @importFrom SummarizedExperiment assay colData assayNames 
#' @importFrom viridis scale_fill_viridis scale_color_viridis
#' @importFrom utils packageVersion write.table
#' @importMethodsFrom BiocGenerics counts updateObject subset colnames rownames
#' @importMethodsFrom scran computeSumFactors
#' @import ggplot2
#' @import SingleCellExperiment
#' @import KernSmooth
NULL
