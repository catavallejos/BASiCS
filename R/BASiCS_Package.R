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
#' @importFrom grDevices adjustcolor colorRampPalette nclass.FD
#' @importFrom ggExtra ggMarginal
#' @importFrom graphics segments
#' @importFrom MASS mvrnorm
#' @importFrom matrixStats colMeans2 colMedians colVars
#' @importFrom matrixStats rowMeans2 rowMedians rowVars   
#' @importFrom methods .hasSlot is new show slotNames Summary 
#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom S4Vectors DataFrame
#' @importFrom scran calculateSumFactors 
#' @importFrom stats acf median model.matrix rgamma rpois runif var lm quantile
#'                   rlnorm
#' @importFrom SingleCellExperiment SingleCellExperiment counts 
#' @importFrom SingleCellExperiment altExp altExp<- 
#' @importFrom SingleCellExperiment altExpNames altExpNames<-
#' @importFrom stats4 plot
#' @importFrom SummarizedExperiment SummarizedExperiment assay assays assayNames
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<-
#' @importFrom viridis scale_color_viridis scale_fill_viridis 
#' @importFrom hexbin hexbin
#' @importFrom utils packageVersion read.delim write.table write.table
#' @importFrom assertthat assert_that
#' @importMethodsFrom BiocGenerics counts colnames rownames subset updateObject 
#' @importMethodsFrom Matrix t rowMeans
#' @importClassesFrom Biobase Versioned
NULL

