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
      ggplot2::labs(title = "Offset MCMC chain", y = "Offset estimate", x = NULL)
    if (Print) {
      print(g)
    }
    Plots[[Type]] <- g
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
    Plots[[Type]] <- g
  }



  if ("MAPlot" %in% Type) {

    # MA plot pre/after offset
    x <- log2(OffsetCorrected@MuBase_old)
    y <- OffsetCorrected@MedianTau_old
    g1 <- ggplot2::ggplot(mapping = ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_hex(bins = NClassFD2D(x, y), 
        aes(fill = ..density..)
      ) + 
      # ggplot2::geom_point() +
      ggplot2::labs(
        x = "Mean expresssion (log2)",
        y = paste(
          "Log2 fold change", OffsetCorrected@GroupLabel1,
          "vs", OffsetCorrected@GroupLabel2),
        title = "Before correction") +
      viridis::scale_fill_viridis(name = "Density") +
      ggplot2::geom_hline(yintercept = 0, lty = 2) +
      ggplot2::geom_hline(yintercept = log2(OffsetCorrected@OffsetEst), lty = 1, col = "red")


    x <- log2(OffsetCorrected@MuBase)
    y <- OffsetCorrected@MedianTau
    g2 <- ggplot2::ggplot(mapping = ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_hex(bins = NClassFD2D(x, y), aes(fill = ..density..)) + 
      # ggplot2::geom_point() +
      ggplot2::labs(
        x = "Mean expresssion (log2)",
        y = paste(
          "Log2 fold change", OffsetCorrected@GroupLabel1,
          "vs", OffsetCorrected@GroupLabel2),
        title = "After correction") +
      ggplot2::geom_hline(yintercept = 0, lty = 2) +
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

NClassFD2D <- function(x, y) {
  max(nclass.FD(x), nclass.FD(y))
}


DEPlots_ResultDE <- function(
    ResultDE,
    GroupLabel1 = ResultDE@GroupLabel1,
    GroupLabel2 = ResultDE@GroupLabel2,
    ProbThresholds = seq(0.5, 0.9995, by = 0.00025),
    Epsilon = ResultDE@Epsilon,
    EFDR = ResultDE@EFDR,
    Table = ResultDE@Table,
    Measure = ResultDE@Name,
    EFDRgrid = ResultDE@EFDRgrid,
    EFNRgrid = ResultDE@EFNRgrid,
    ProbThreshold = ResultDE@ProbThreshold,
    Which = c("MAPlot", "VolcanoPlot", "GridPlot")){

  DEPlots(
    GroupLabel1 = GroupLabel1,
    GroupLabel2 = GroupLabel2,
    ProbThresholds = ProbThresholds,
    Epsilon = Epsilon,
    EFDR = EFDR,
    Table = Table,
    Measure = Measure,
    EFDRgrid = EFDRgrid,
    EFNRgrid = EFNRgrid,
    ProbThreshold = ProbThreshold,
    Which = Which
  )
}

DEPlots <- function(
    GroupLabel1,
    GroupLabel2,
    ProbThresholds = seq(0.5, 0.9995, by = 0.00025),
    Epsilon,
    EFDR,
    Table,
    Measure,
    EFDRgrid,
    EFNRgrid,
    ProbThreshold,
    Which = c("MAPlot", "VolcanoPlot", "GridPlot")
    ) {

  Which <- match.arg(Which, several.ok = TRUE)

  Plots <- list()

  if ("MAPlot" %in% Which) {
    Plots <- c(Plots, 
      list(MAPlot(Measure, Table, GroupLabel1, GroupLabel2, Epsilon))
    )
  }
  if ("VolcanoPlot" %in% Which) {
    Plots <- c(Plots, 
      list(VolcanoPlot(Measure, Table, GroupLabel1, GroupLabel2, Epsilon))
    )
  }
  if ("GridPlot" %in% Which) {
    Plots <- c(
      list(GridPlot(Measure, EFDR, ProbThresholds, EFDRgrid, EFNRgrid, ProbThreshold)),
      Plots
    )
  }
  if (length(Plots) == 1) {
    Plots <- Plots[[1]]
  } else {
    Plots <- cowplot::plot_grid(plotlist = Plots, nrow = 1)
  }
  Plots
}



GridPlot <- function(Measure, 
                       EFDR,
                       ProbThresholds = seq(0.5, 0.9995, by = 0.00025),
                       EFDRgrid,
                       EFNRgrid,
                       ProbThreshold
                       ) {
  df <- data.frame(
    ProbThresholds, 
    EFDR = EFDRgrid, 
    EFNR = EFNRgrid)

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
    ggplot2::geom_vline(xintercept = ProbThreshold, col = "indianred", lty = 1)
}



MAPlot <- function(Measure, Table, GroupLabel1, GroupLabel2, Epsilon) {
  IndDiff <- !Table[[paste0("ResultDiff", Measure)]] %in% c("ExcludedByUser", "NoDiff")
  bins <- NClassFD2D(
    Table[[paste0(Measure, "Overall")]],
    Table[[paste0(Measure, DistanceVar(Measure))]]
  )
  ggplot2::ggplot(
      Table, 
      ggplot2::aes_string(
        x = paste0(Measure, "Overall"), 
        y = paste0(Measure, DistanceVar(Measure)))
    ) + 
    ggplot2::geom_hex(bins = bins, aes(fill = ..density..)) +
    # ggplot2::geom_point() +
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
  bins <- NClassFD2D(
    Table[[paste0(Measure, DistanceVar(Measure))]],
    Table[[paste0("ProbDiff", Measure)]]
  )
  ggplot2::ggplot(
      Table, 
      ggplot2::aes_string(
        x = paste0(Measure, DistanceVar(Measure)), 
        y = paste0("ProbDiff", Measure))
    ) +
    # ggplot2::geom_point() +
    ggplot2::geom_hex(bins = bins, aes(fill = ..density..)) +
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


BASiCS_PlotDE <- function(ResultsDE, Which = c("MAPlot", "VolcanoPlot", "GridPlot"), ...) {
  l <- lapply(ResultsDE@Results, DEPlots_ResultDE, Which = Which, ...)
  nrow <- if (length(Which) > 1) length(ResultsDE@Results) else length(Which)
  cowplot::plot_grid(plotlist = l, nrow = nrow)
}
