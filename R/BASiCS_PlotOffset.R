#' @name BASiCS_PlotOffset
#'
#' @title Visualise global offset in mean expression between two chains.
#'
#' @description Visualise global offset in mean expression between two
#' \code{BASiCS_Chain} objects.
#'
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
#' data("ChainSC")
#' data("ChainRNA")
#'
#' Corrected <- BASiCS_CorrectOffset(ChainSC, ChainRNA, "a", "b", Plot = FALSE)
#' BASiCS_PlotOffset(Corrected)
#'
#' @return Plot objects.
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#' @author Alan O'Callaghan \email{a.b.o'callaghan@sms.ed.ac.uk}
#'
#' @export
BASiCS_PlotOffset <- function(Chain1,
                              Chain2,
                              GroupLabel1 = "Group 1",
                              GroupLabel2 = "Group 2",
                              Type = c(
                                "offset estimate",
                                "before-after",
                                "MAPlot"
                              ),
                              Print = FALSE
                              ) {
  L <- BASiCS_CorrectOffset(Chain1, Chain2)

  OffsetChain <- L$OffsetChain
  OffsetEst <- L$OffsetEst
  Chain1_offset <- L$Chain
  Chain1_offset@parameters$mu <- Chain1@parameters$mu / OffsetEst
  Chain2_offset <- Chain2  # Chain2 requires no change

  Mu1 <- matrixStats::colMedians(Chain1_offset@parameters$mu)
  Mu2 <- matrixStats::colMedians(Chain2@parameters$mu)
  Mu1_old <- matrixStats::colMedians(Chain1@parameters$mu)

  Type <- match.arg(Type)

  Plots <- list()
  if ("offset estimate" %in% Type) {
    # Offset uncertainty
    g <- ggplot2::ggplot(mapping = ggplot2::aes(y = OffsetChain)) +
      ggplot2::geom_boxplot(na.rm = TRUE) +
      ggplot2::labs(
        # title = "Offset MCMC chain",
        y = "Offset estimate", x = NULL) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
    if (Print) {
      print(g)
    }
    Plots[[Type]] <- g
  }

  # Mean expression parameters before/after offset correction
  if ("before-after" %in% Type) {
    ylim <- range(
      c(Mu1_old, Mu1, Mu2)
    )
    df <- cbind(Mu1_old, Mu2)
    colnames(df) <- c(GroupLabel1, GroupLabel2)
    mdf <- reshape2::melt(df)
    g1 <- ggplot2::ggplot(mdf, ggplot2::aes_string(x = "Var2", y = "value")) +
      ggplot2::labs(
        title = "Before correction",
        y = "Mean expression",
        x = NULL
      ) +
      ggplot2::geom_violin(na.rm = TRUE) +
      ggplot2::scale_y_log10(limits = ylim)

    df <- cbind(Mu1, Mu2)
    colnames(df) <- c(GroupLabel1, GroupLabel2)
    mdf <- reshape2::melt(df)
    g2 <- ggplot2::ggplot(mdf, ggplot2::aes_string(x = "Var2", y = "value")) +
      ggplot2::labs(
        title = "After correction",
        y = "Mean expression",
        x = NULL
      ) +
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
      x = log2(MuBase_old),
      y = MedianTau_old
    )
    g1 <- ggplot2::ggplot(df, mapping = ggplot2::aes_string(x = "x", y = "y")) +
      ggplot2::geom_hex(
        # bins = NClassFD2D(df$x, df$y),
        bins = 100,
        aes_string(fill = "..density.."),
        na.rm = TRUE) +
      # ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0, lty = 2) +
      ggplot2::geom_hline(
        yintercept = log2(OffsetEst),
        lty = 1,
        col = "red",
        na.rm = TRUE) +
      ggplot2::labs(
        x = "Mean expresssion (log2)",
        y = paste(
          "Log2 fold change", GroupLabel1,
          "vs", GroupLabel2),
        title = "Before correction") +
      viridis::scale_fill_viridis(name = "Density", guide = FALSE)

    df <- data.frame(
      x = log2(MuBase),
      y = MedianTau
    )
    g2 <- ggplot2::ggplot(df, mapping = ggplot2::aes_string(x = "x", y = "y")) +
      ggplot2::geom_hex(
        mapping = aes_string(fill = "..density.."),
        bins = 100,
        # bins = NClassFD2D(df$x, df$y),
        bins = 100,
        na.rm = TRUE) +
      # ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0, lty = 2) +
      ggplot2::labs(
        x = "Mean expresssion (log2)",
        y = paste(
          "Log2 fold change", GroupLabel1,
          "vs", GroupLabel2),
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
  } else {
    Plots <- cowplot::plot_grid(plotlist = Plots)
  }
}
