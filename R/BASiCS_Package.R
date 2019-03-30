#' @useDynLib BASiCS
#' @exportPattern "^[^\\Hidden]"
#' @importFrom coda HPDinterval mcmc effectiveSize
#' @importFrom cowplot plot_grid
#' @importFrom data.table fread
#' @importFrom ggplot2 aes aes_string 
#' @importFrom ggplot2 geom_ribbon geom_segment ggplot ggtitle 
#'                     geom_hline geom_line geom_point 
#' @importFrom ggplot2 scale_fill_gradientn scale_x_log10 scale_y_log10 
#'                     scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 theme_minimal theme_classic 
#' @importFrom ggplot2 xlab ylab labs
#' @importFrom graphics abline barplot legend lines par plot points segments 
#' @importFrom graphics smoothScatter 
#' @importFrom grDevices adjustcolor colorRampPalette
#' @importFrom ggExtra ggMarginal
#' @importFrom MASS mvrnorm
#' @importFrom matrixStats colMeans2 colMedians colVars
#' @importFrom matrixStats rowMeans2 rowMedians rowVars   
#' @importFrom methods .hasSlot is new show slotNames Summary 
#' @importFrom Rcpp evalCpp
#' @importFrom S4Vectors DataFrame metadata 
#' @importFrom stats acf median model.matrix rgamma rpois runif 
#' @importFrom SingleCellExperiment isSpike isSpike<- counts clearSpikes
#' @importFrom SummarizedExperiment assay assayNames  colData 
#' @importFrom viridis scale_color_viridis scale_fill_viridis 
#' @importFrom utils packageVersion write.table
#' @importMethodsFrom BiocGenerics counts colnames rownames subset updateObject 
#' @importMethodsFrom scran computeSumFactors
#' @importClassesFrom Biobase Versioned
NULL
