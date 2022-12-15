#' @name BASiCS_PlotOffset
#'
#' @title Visualise global offset in mean expression between two chains.
#'
#' @description Visualise global offset in mean expression between two
#' \code{BASiCS_Chain} objects.
#'
#' @param Chain1,Chain2 \linkS4class{BASiCS_Chain} objects to be plotted.
#' @param Type The type of plot generated.
#'    \code{"offset estimate"} produces a boxplot of the offset alongside
#'    an estimate of the global offset.
#'    \code{"before-after"} produces MA plots of Mean expression against
#'    log2(fold-change) before and after offset correction.
#'    \code{"MA plot"} produces an MA plot of Mean expression against
#'    log2(fold-change).
#' @param GroupLabel1,GroupLabel2 Labels for Chain1 and Chain2 in the resulting
#' plot(s).
#'
#' @examples
#'
#' # Loading two 'BASiCS_Chain' objects (obtained using 'BASiCS_MCMC')
#' data("ChainSC")
#' data("ChainRNA")
#'
#' BASiCS_PlotOffset(ChainSC, ChainRNA)
#'
#' @return Plot objects.
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#' @author Alan O'Callaghan
#'
#' @export
BASiCS_PlotOffset <- function(Chain1,
                              Chain2,
                              Type = c(
                                "offset estimate",
                                "before-after",
                                "MAPlot"
                              ),
                              GroupLabel1 = "Group 1",
                              GroupLabel2 = "Group 2"
                              ) {
  L <- BASiCS_CorrectOffset(Chain1, Chain2)
  
  OffsetChain <- L$OffsetChain
  OffsetEst <- L$Offset
  Chain1_offset <- L$Chain
  Chain2_offset <- Chain2  # Chain2 requires no change

  Mu1 <- matrixStats::colMedians(Chain1_offset@parameters$mu)
  Mu2 <- matrixStats::colMedians(Chain2@parameters$mu)
  Mu1_old <- matrixStats::colMedians(Chain1@parameters$mu)
  
  n1 <- ncol(Chain1@parameters$nu)
  n2 <- ncol(Chain2@parameters$nu)
  n <- n1 + n2
  MuBase <- (Mu1 * n1 + Mu2 * n2) / n

  Type <- match.arg(Type, several.ok = TRUE)

  Plots <- list()
  if ("offset estimate" %in% Type) {
    # Offset uncertainty
    g <- ggplot(mapping = aes(y = OffsetChain)) +
      geom_boxplot(na.rm = TRUE) +
      labs(
        # title = "Offset MCMC chain",
        y = "Offset estimate", x = NULL) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    Plots <- c(Plots, list(g))
  }

  # Mean expression parameters before/after offset correction
  if ("before-after" %in% Type) {
    ylim <- range(
      c(Mu1_old, Mu1, Mu2)
    )
    df <- cbind(Mu1_old, Mu2)
    colnames(df) <- c(GroupLabel1, GroupLabel2)
    mdf <- reshape2::melt(df)
    g1 <- ggplot(mdf, aes(x = .data$Var2, y = .data$value)) +
      labs(
        title = "Before correction",
        y = "Mean expression",
        x = NULL
      ) +
      geom_violin(na.rm = TRUE) +
      scale_y_log10(limits = ylim)

    df <- cbind(Mu1, Mu2)
    colnames(df) <- c(GroupLabel1, GroupLabel2)
    mdf <- reshape2::melt(df)
    g2 <- ggplot(mdf, aes(x = .data$Var2, y = .data$value)) +
      labs(
        title = "After correction",
        y = "Mean expression",
        x = NULL
      ) +
      geom_violin() +
      scale_y_log10(limits = ylim)

    g <- cowplot::plot_grid(g1, g2)
    Plots <- c(Plots, list(g))
  }

  if ("MAPlot" %in% Type) {
    MuBase_old <- (Mu1_old * n1 + Mu2 * n2) / n
    ChainTau_old <- log2(Chain1@parameters$mu / Chain2@parameters$mu)
    MedianTau_old <- matrixStats::colMedians(ChainTau_old)
    ChainTau <- log2(Chain1_offset@parameters$mu / Chain2_offset@parameters$mu)
    MedianTau <- matrixStats::colMedians(ChainTau)

    # MA plot pre/after offset
    df <- data.frame(
      x = log2(MuBase_old),
      y = MedianTau_old
    )
    g1 <- ggplot(df, mapping = aes(x = .data$x, y = .data$y)) +
      geom_hex(
        # bins = NClassFD2D(df$x, df$y),
        bins = 100,
        aes(fill = after_stat(density)),
        na.rm = TRUE) +
      # geom_point() +
      geom_hline(yintercept = 0, lty = 2) +
      geom_hline(
        yintercept = log2(OffsetEst),
        lty = 1,
        col = "red",
        na.rm = TRUE) +
      labs(
        x = "Mean expresssion (log2)",
        y = paste(
          "Log2 fold change", GroupLabel1,
          "vs", GroupLabel2),
        title = "Before correction") +
      scale_fill_viridis(name = "Density", guide = "none")

    df <- data.frame(
      x = log2(MuBase),
      y = MedianTau
    )
    g2 <- ggplot(df, mapping = aes(x = .data$x, y = .data$y)) +
      geom_hex(
        mapping = aes(fill = after_stat(density)),
        bins = 100,
        # bins = NClassFD2D(df$x, df$y),
        na.rm = TRUE) +
      # geom_point() +
      geom_hline(yintercept = 0, lty = 2) +
      labs(
        x = "Mean expresssion (log2)",
        y = paste(
          "Log2 fold change", GroupLabel1,
          "vs", GroupLabel2),
        title = "After correction") +
      scale_fill_viridis(name = "Density")

    g <- cowplot::plot_grid(g1, g2)
    Plots <- c(Plots, list(g))
  }
  if (length(Plots) == 1) {
    Plots <- Plots[[1]]
  } else {
    Plots <- cowplot::plot_grid(plotlist = Plots, nrow = 1)
  }
  Plots
}
