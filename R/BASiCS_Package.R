#' @useDynLib BASiCS
#' @importFrom coda HPDinterval mcmc effectiveSize
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 aes aes_string 
#' @importFrom ggplot2 geom_ribbon geom_segment ggplot ggtitle 
#'                     geom_hline geom_line geom_point 
#' @importFrom ggplot2 scale_fill_gradientn scale_x_log10 scale_y_log10 
#'                     scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 theme_minimal theme_classic 
#' @importFrom ggplot2 xlab ylab labs
#' @importFrom graphics abline barplot legend lines par points segments 
#' @importFrom graphics smoothScatter 
#' @importFrom grDevices adjustcolor colorRampPalette nclass.FD
#' @importFrom ggExtra ggMarginal
#' @importFrom KernSmooth bkde2D
#' @importFrom MASS mvrnorm
#' @importFrom matrixStats colMeans2 colMedians
#' @importFrom matrixStats rowMeans2 rowMedians
#' @importFrom Matrix t rowSums colSums rowMeans colMeans 
#' @importFrom methods .hasSlot is new show slotNames Summary 
#' @importFrom Rcpp evalCpp
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' @importFrom scran calculateSumFactors
#' @importFrom stats acf median model.matrix rgamma rpois runif var
#' @importFrom SingleCellExperiment SingleCellExperiment counts 
#' @importFrom SingleCellExperiment altExp altExp<- 
#' @importFrom SingleCellExperiment altExpNames altExpNames<-
#' @importFrom stats4 plot
#' @importFrom SummarizedExperiment SummarizedExperiment assay assayNames
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom viridis scale_color_viridis scale_fill_viridis 
#' @importFrom hexbin hexbin
#' @importFrom utils packageVersion read.delim write.table
#' @importMethodsFrom BiocGenerics counts colnames rownames subset updateObject 
#' @importClassesFrom Biobase Versioned
NULL

