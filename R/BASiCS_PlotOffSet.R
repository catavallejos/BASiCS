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
                                "MA plot"
                              ),
                              Print = FALSE
                              ) {
  Type <- match.arg(Type)

  plots <- list()
  if ("offset estimate" %in% Type) {
    # Offset uncertainty
    g <- ggplot2::ggplot(mapping = ggplot2::aes(y = OffsetCorrected@OffsetChain)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(title = "Offset MCMC chain", y = "Offset estimate", x = NULL)
    if (Print) {
      print(g)
    }
    plots[[Type]] <- g
  }

  # Mean expression parameters before/after offset correction
  if ("before-after" %in% Type) {
    df <- cbind(OffsetCorrected@Mu1_old, OffsetCorrected@Mu2)
    colnames(df) <- c(OffsetCorrected@GroupLabel1, OffsetCorrected@GroupLabel2)
    mdf <- reshape2::melt(df)
    ggplot2::ggplot(mdf, ggplot2::aes(x = Var2, y = value)) +
      ggplot2::labs(title = "Before correction", y = "Mean expression", x = NULL) +
      ggplot2::geom_boxplot() +
      ggplot2::scale_y_log10()

    df <- cbind(OffsetCorrected@Mu1, OffsetCorrected@Mu2)
    colnames(df) <- c(OffsetCorrected@GroupLabel1, OffsetCorrected@GroupLabel2)
    mdf <- reshape2::melt(df)
    g <- ggplot2::ggplot(mdf, ggplot2::aes(x = Var2, y = value)) +
      ggplot2::labs(title = "After correction", y = "Mean expression", x = NULL) +
      ggplot2::geom_boxplot() +
      ggplot2::scale_y_log10()

    if (Print) {
      print(g)
    }
    plots[[Type]] <- g
  }

  if ("MA plot" %in% Type) {
    # MA plot pre/after offset

    g1 <- ggplot2::ggplot(mapping = ggplot2::aes(x = log2(OffsetCorrected@MuBase_old), 
                                                 y = OffsetCorrected@MedianTau_old)) +
      # geom_hex() + 
      ggplot2::geom_point() +
      ggplot2::labs(
        x = "Mean expresssion (log2)",
        y = paste(
          "Log2 fold change", OffsetCorrected@GroupLabel1,
          "vs", OffsetCorrected@GroupLabel2),
        title = "Before correction") +
      viridis::scale_fill_viridis() +
      ggplot2::geom_hline(yintercept = 0, lty = 2) +
      ggplot2::geom_hline(yintercept = log2(OffsetCorrected@OffsetEst), lty = 1, col = "red")

    g2 <- ggplot2::ggplot(mapping = ggplot2::aes(x = log2(OffsetCorrected@MuBase), 
                                                 y = OffsetCorrected@MedianTau)) +
      # ggplot2::geom_hex() + 
      ggplot2::geom_point() +
      ggplot2::labs(
        x = "Mean expresssion (log2)",
        y = paste(
          "Log2 fold change", OffsetCorrected@GroupLabel1,
          "vs", OffsetCorrected@GroupLabel2),
        title = "After correction") +
      ggplot2::geom_hline(yintercept = 0, lty = 2) +
      viridis::scale_fill_viridis()

    g <- cowplot::plot_grid(g1, g2)
    if (Print) {
      print(g)
    }
    plots[[Type]] <- g
  }
  plots
}




DEPlots <- function(
    GroupLabel1,
    GroupLabel2,
    ProbThresholds = seq(0.5, 0.9995, by = 0.00025),
    Aux,
    Epsilon,
    EFDR,
    Table,
    Measure,
    Search = TRUE
    ) {

  plots <- list(
    MAPlot(Measure, Table, GroupLabel1, GroupLabel2, Epsilon),
    VolcanoPlot(Measure, Table, GroupLabel1, GroupLabel2, Epsilon)
  )

  if (Search) {
    plots <- c(
      list(PlotSearch(Measure, Aux, EFDR, ProbThresholds)),
      plots
    )
  }

  cowplot::plot_grid(plotlist = plots, nrow = 1)
}



PlotSearch <- function(Measure, 
                       Aux, 
                       EFDR,
                       ProbThresholds = seq(0.5, 0.9995, by = 0.00025)
                       ) {
  df <- data.frame(
    ProbThresholds, 
    EFDR = Aux$EFDRgrid, 
    EFNR = Aux$EFNRgrid)

  mdf <- reshape2::melt(df,measure.vars = c("EFDR", "EFNR"))
  ggplot2::ggplot(mdf, aes(x = ProbThresholds, y = value, color = variable)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      y = "Error rate", 
      x = "Probability threshold",
      title = paste("Differential", MeasureName(Measure))
    ) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::scale_color_discrete(name = ""
      # , palette = "Dark2"
      ) +
    ggplot2::geom_hline(yintercept = EFDR, col = "steelblue", lty = 2) +
    ggplot2::geom_vline(xintercept = Aux$OptThreshold[[1]], col = "indianred", lty = 1)
}



PlotGrids <- function(
    GroupLabel1,
    GroupLabel2,
    ProbThresholds = seq(0.5, 0.9995, by = 0.00025), 
    AuxMean, 
    AuxDisp, 
    AuxResDisp = NULL,
    EpsilonM,
    EpsilonD,
    EpsilonR,
    EFDR_M,
    EFDR_D,
    EFDR_R = NULL,
    TableMean,
    TableDisp,
    TableResDisp = NULL,
    IncludeEpsilon = !is.null(AuxResDisp)
    ) {

  plots <- list(
    PlotSearch("Mean", AuxMean, EFDR_M, ProbThresholds),
    PlotSearch("Disp", AuxDisp, EFDR_D, ProbThresholds)
  )

  if (IncludeEpsilon) {
    plots <- c(
      plots, 
      list(PlotSearch("ResDisp", AuxResDisp, EFDR_R, ProbThresholds))
    )
  }
  cowplot::plot_grid(plotlist = plots, nrow = 1)
}



PlotMAs <- function(
    GroupLabel1,
    GroupLabel2,
    ProbThresholds = seq(0.5, 0.9995, by = 0.00025), 
    AuxMean, 
    AuxDisp, 
    AuxResDisp = NULL,
    EpsilonM,
    EpsilonD,
    EpsilonR,
    EFDR_M,
    EFDR_D,
    EFDR_R = NULL,
    TableMean,
    TableDisp,
    TableResDisp = NULL,
    IncludeEpsilon = !is.null(AuxResDisp)
    ) {

  # MA plots
  plots <- list(
    MAPlot("Mean", TableMean, GroupLabel1, GroupLabel2, EpsilonM),
    MAPlot("Disp", TableDisp, GroupLabel1, GroupLabel2, EpsilonD)
  )

  if (IncludeEpsilon) {
    plots <- c(
      plots, 
      list(MAPlot("ResDisp", TableResDisp, GroupLabel1, GroupLabel2, EpsilonR))
    )
  }
  cowplot::plot_grid(plotlist = plots, nrow = 1)
}



PlotVolcanos <- function(
    GroupLabel1,
    GroupLabel2,
    ProbThresholds = seq(0.5, 0.9995, by = 0.00025), 
    AuxMean, 
    AuxDisp, 
    AuxResDisp = NULL,
    EpsilonM,
    EpsilonD,
    EpsilonR,
    EFDR_M,
    EFDR_D,
    EFDR_R = NULL,
    TableMean,
    TableDisp,
    TableResDisp = NULL,
    IncludeEpsilon = !is.null(AuxResDisp)
    ) {

  # MA plots
  plots <- list(
    VolcanoPlot("Mean", TableMean, GroupLabel1, GroupLabel2, EpsilonM),
    VolcanoPlot("Disp", TableDisp, GroupLabel1, GroupLabel2, EpsilonD)
  )

  if (IncludeEpsilon) {
    plots <- c(
      plots, 
      list(VolcanoPlot("ResDisp", TableResDisp, GroupLabel1, GroupLabel2, EpsilonR))
    )
  }
  cowplot::plot_grid(plotlist = plots, nrow = 1)

}

MAPlot <- function(Measure, Table, GroupLabel1, GroupLabel2, Epsilon) {
  IndDiff <- !Table[[paste0("ResultDiff", Measure)]] %in% c("ExcludedByUser", "NoDiff")
  ggplot2::ggplot(
      Table, 
      ggplot2::aes_string(
        x = paste0(Measure, "Overall"), 
        y = paste0(Measure, DistanceVar(Measure)))
    ) + 
    ggplot2::geom_hex(bins = 75, aes(fill = ..density..)) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(trans = "log2") +
    viridis::scale_fill_viridis(name = "Density") +
    ggplot2::labs(
      x = paste(cap(MeasureName(Measure)), " (log2)"),
      ylab = paste("Log2 fold change",
        GroupLabel1, "vs",
        GroupLabel2),
      title = paste("Differential", MeasureName(Measure))
    ) +
    ggplot2::geom_point(data = Table[IndDiff, ], shape = 16, colour = "violetred") +
    ggplot2::geom_hline(yintercept = -Epsilon, lty = "dashed", color = "grey40") +
    ggplot2::geom_hline(yintercept = Epsilon, lty = "dashed", color = "grey40")
}

VolcanoPlot <- function(Measure, Table, GroupLabel1, GroupLabel2, Epsilon) {
  IndDiff <- !Table[[paste0("ResultDiff", Measure)]] %in% c("ExcludedByUser", "NoDiff")
  ggplot2::ggplot(
      Table, 
      ggplot2::aes_string(
        x = paste0(Measure, DistanceVar(Measure)), 
        y = paste0("ProbDiff", Measure))
    ) +
    # ggplot2::geom_point() +
    ggplot2::geom_hex(bins = 75, aes(fill = ..density..)) +
    viridis::scale_fill_viridis(name = "Density") +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::labs(y = "Posterior probability",
         x = paste("Log2 fold change", GroupLabel1, "vs", GroupLabel2),
         title = paste("Differential", MeasureName(Measure), "test")
    ) +
    ggplot2::geom_point(data = Table[IndDiff, ], shape = 16, col = "violetred") +
    ggplot2::geom_hline(yintercept = -Epsilon, lty = "dashed", color = "grey40") +
    ggplot2::geom_hline(yintercept = Epsilon, lty = "dashed", color = "grey40")
}


BASiCS_PlotDE <- function(
    GroupLabel1,
    GroupLabel2,
    ProbThresholds = seq(0.5, 0.9995, by = 0.00025), 
    AuxMean, 
    AuxDisp, 
    AuxResDisp = NULL,
    EpsilonM,
    EpsilonD,
    EpsilonR,
    EFDR_M,
    EFDR_D,
    EFDR_R = NULL,
    TableMean,
    TableDisp,
    TableResDisp = NULL,
    IncludeEpsilon = !is.null(AuxResDisp),
    Search = TRUE) {

  if (Search) {

    ProbThresholds <- seq(0.5, 0.9995, by = 0.00025)

    plots <- list(
      PlotSearch("Mean", AuxMean, EFDR_M, ProbThresholds),
      PlotSearch("Disp", AuxDisp, EFDR_D, ProbThresholds)
    )

    if (IncludeEpsilon) {
      plots <- c(
        plots, 
        list(PlotSearch("ResDisp", AuxResDisp, EFDR_R, ProbThresholds))
      )
    }
    c <- cowplot::plot_grid(plotlist = plots, nrow = 1)
    print(c)
  }

  # MA plots
  plots <- list(
    MAPlot("Mean", TableMean, GroupLabel1, GroupLabel2, EpsilonM),
    MAPlot("Disp", TableDisp, GroupLabel1, GroupLabel2, EpsilonD)
  )

  if (IncludeEpsilon) {
    plots <- c(plots, list(MAPlot("ResDisp", TableResDisp, GroupLabel1, GroupLabel2, EpsilonR)))
  }
  c <- cowplot::plot_grid(plotlist = plots, nrow = 1)
  print(c)

  plots <- list(
    VolcanoPlot("Mean", TableMean, GroupLabel1, GroupLabel2, EpsilonM),
    VolcanoPlot("Disp", TableDisp, GroupLabel1, GroupLabel2, EpsilonD)
  )

  if (IncludeEpsilon) {
    plots <- c(
      plots, 
      list(VolcanoPlot("ResDisp", TableResDisp, GroupLabel1, GroupLabel2, EpsilonR))
    )
  }
  c <- cowplot::plot_grid(plotlist = plots, nrow = 1)
  print(c)
}
