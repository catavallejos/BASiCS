#' @name BASiCS_PlotOffset
#'
#' @title Visualise global offset in mean expression between two chains.
#'
#' @description Visualise global offset in mean expression between two 
#' \code{BASiCS_Chain} objects.
#' 
#' @param OffsetCorrected Object of class BASiCS_OffsetCorrected
#' @param Type The type of plot generated.
#'    \code{"offset estimate"} produces a boxplot of the offset alongside
#'    an estimate of the global offset.
#'    \code{"before-after"} produces MA plots of Mean expression against
#'    log2(fold-change) before and after offset correction.
#'    \code{"MA plot"} produces an MA plot of Mean expression against 
#'    log2(fold-change).
#' @param Print Print the plots (TRUE), or simply return them as a list (FALSE)?
#' Default is \code{FALSE}.
#'
#' @examples
#'
#' # Loading two 'BASiCS_Chain' objects (obtained using 'BASiCS_MCMC')
#' data(ChainSC)
#' data(ChainRNA)
#' 
#' Corrected <- BASiCS_CorrectOffset(ChainSC, ChainRNA, "a", "b", Plot = FALSE)
#' BASiCS_PlotOffset(Corrected)
#'
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#' @author Alan O'Callaghan \email{a.b.o'callaghan@sms.ed.ac.uk}
#' 
#' @export
BASiCS_PlotOffset <- function(OffsetCorrected,
                              Type = c(
                                "offset estimate", 
                                "before-after",
                                "MAPlot"
                              ),
                              Print = FALSE
                              ) {
  Type <- match.arg(Type)

  Plots <- list()
  if ("offset estimate" %in% Type) {
    # Offset uncertainty
    g <- ggplot2::ggplot(mapping = ggplot2::aes(y = OffsetCorrected@OffsetChain)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(
        # title = "Offset MCMC chain", 
        y = "Offset estimate", x = NULL) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    if (Print) {
      print(g)
    }
    Plots[[Type]] <- g
  }

  # Mean expression parameters before/after offset correction
  if ("before-after" %in% Type) {
    ylim <- range(
      c(OffsetCorrected@Mu1_old, OffsetCorrected@Mu1, OffsetCorrected@Mu2)
    )
    df <- cbind(OffsetCorrected@Mu1_old, OffsetCorrected@Mu2)
    colnames(df) <- c(OffsetCorrected@GroupLabel1, OffsetCorrected@GroupLabel2)
    mdf <- reshape2::melt(df)
    g1 <- ggplot2::ggplot(mdf, ggplot2::aes(x = Var2, y = value)) +
      ggplot2::labs(title = "Before correction", y = "Mean expression", x = NULL) +
      ggplot2::geom_violin() +
      ggplot2::scale_y_log10(limits = ylim)

    df <- cbind(OffsetCorrected@Mu1, OffsetCorrected@Mu2)
    colnames(df) <- c(OffsetCorrected@GroupLabel1, OffsetCorrected@GroupLabel2)
    mdf <- reshape2::melt(df)
    g2 <- ggplot2::ggplot(mdf, ggplot2::aes(x = Var2, y = value)) +
      ggplot2::labs(title = "After correction", y = "Mean expression", x = NULL) +
      ggplot2::geom_violin() +
      ggplot2::scale_y_log10(limits = ylim)

    g <- cowplot::plot_grid(g1, g2)
    if (Print) {
      print(g)
    }
    Plots[[Type]] <- g
  }



  if ("MAPlot" %in% Type) {

    # MA plot pre/after offset
    df <- data.frame(
      x = log2(OffsetCorrected@MuBase_old),
      y = OffsetCorrected@MedianTau_old
    )
    g1 <- ggplot2::ggplot(df, mapping = ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_hex(
        # bins = NClassFD2D(df$x, df$y), 
        bins = 100,
        aes(fill = ..density..), 
        na.rm = TRUE) + 
      # ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0, lty = 2) +
      ggplot2::geom_hline(
        yintercept = log2(OffsetCorrected@OffsetEst), 
        lty = 1, 
        col = "red", 
        na.rm = TRUE) +
      ggplot2::labs(
        x = "Mean expresssion (log2)",
        y = paste(
          "Log2 fold change", OffsetCorrected@GroupLabel1,
          "vs", OffsetCorrected@GroupLabel2),
        title = "Before correction") +
      viridis::scale_fill_viridis(name = "Density")

    df <- data.frame(
      x = log2(OffsetCorrected@MuBase),
      y = OffsetCorrected@MedianTau
    )
    g2 <- ggplot2::ggplot(df, mapping = ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_hex(
        mapping = aes(fill = ..density..), 
        bins = 100,
        # bins = NClassFD2D(df$x, df$y), 
        na.rm = TRUE) + 
      # ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0, lty = 2) +
      ggplot2::labs(
        x = "Mean expresssion (log2)",
        y = paste(
          "Log2 fold change", OffsetCorrected@GroupLabel1,
          "vs", OffsetCorrected@GroupLabel2),
        title = "After correction") +
      viridis::scale_fill_viridis(name = "Density")

    g <- cowplot::plot_grid(g1, g2)
    if (Print) {
      print(g)
    }
    Plots[[Type]] <- g
  }
  if (length(Plots) == 1) {
    Plots <- Plots[[1]]
  }
  Plots
}

MAPlot <- function(Measure, Table, GroupLabel1, GroupLabel2, Epsilon) {

  IndDiff <- !Table[[paste0("ResultDiff", Measure)]] %in% c("ExcludedByUser", "NoDiff")
  # bins <- NClassFD2D(
  #   Table[[paste0(Measure, "Overall")]],
  #   Table[[paste0(Measure, DistanceVar(Measure))]]
  # )
  bins <- 100
  ggplot2::ggplot(
      Table, 
      ggplot2::aes_string(
        x = paste0(Measure, "Overall"), 
        y = paste0(Measure, DistanceVar(Measure)))
    ) + 
    ggplot2::geom_hex(bins = bins, aes(fill = ..density..), na.rm = TRUE) +
    # ggplot2::geom_point() +
    ggplot2::geom_point(
      data = Table[IndDiff, ], 
      shape = 16, 
      colour = "violetred", 
      na.rm = TRUE) +
    ggplot2::geom_hline(
      yintercept = c(-Epsilon, Epsilon), 
      lty = "dashed", 
      color = "grey40", 
      na.rm = TRUE) +
    ggplot2::scale_x_continuous(trans = "log2") +
    viridis::scale_fill_viridis(name = "Density") +
    ggplot2::labs(
      x = paste(cap(MeasureName(Measure))),
      y = paste(cap(LogDistanceName(Measure)),
        GroupLabel1, "vs",
        GroupLabel2)
      # ,
      # title = paste("Differential", MeasureName(Measure))
    )
}
