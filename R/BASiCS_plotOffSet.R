#' @export
BASiCS_plotOffset <- function(Chain1, 
                              Chain2,
                              GroupLabel1 = "Group1",
                              GroupLabel2 = "Group2",
                              Type = c(
                                "offset estimate", 
                                "before-after",
                                "MA plot"
                              ),
                              Print = FALSE
                              ) {
  Type <- match.arg(Type)

  n1 <- ncol(Chain1@parameters$nu)
  n2 <- ncol(Chain2@parameters$nu)
  n <- n1 + n2

  # Calculating iteration-specific offset
  OffsetChain <- matrixStats::rowSums2(Chain1@parameters$mu) /
                  matrixStats::rowSums2(Chain2@parameters$mu)
  # Offset point estimate
  OffsetEst <- median(OffsetChain)

  # Offset correction
  Chain1_offset <- Chain1
  Chain1_offset@parameters$mu <- Chain1@parameters$mu / OffsetEst
  Chain2_offset <- Chain2  # Chain2 requires no change

  Mu1 <- matrixStats::colMedians(Chain1_offset@parameters$mu)
  Mu2 <- matrixStats::colMedians(Chain2_offset@parameters$mu)

  Delta1 <- matrixStats::colMedians(Chain1_offset@parameters$delta)
  Delta2 <- matrixStats::colMedians(Chain2_offset@parameters$delta)

  Mu1_old <- matrixStats::colMedians(Chain1@parameters$mu)

  MuBase_old <- (Mu1_old * n1 + Mu2 * n2) / n

  ChainTau_old <- log2(Chain1@parameters$mu / Chain2@parameters$mu)
  MedianTau_old <- matrixStats::colMedians(ChainTau_old)

  # Offset corrected LFC estimates
  MuBase <- (Mu1 * n1 + Mu2 * n2) / n
  ChainTau <- log2(Chain1_offset@parameters$mu / Chain2_offset@parameters$mu)
  MedianTau <- matrixStats::colMedians(ChainTau)

  plots <- list()
  if ("offset estimate" %in% Type) {
    # Offset uncertainty
    g <- ggplot2::ggplot(mapping = aes(y = OffsetChain)) +
      geom_boxplot() +
      labs(title = "Offset MCMC chain", y = "Offset estimate", x = NULL)
    if (Print) {
      print(g)
    }
    plots[[Type]] <- g
  }

  # Mean expression parameters before/after offset correction
  if ("before-after" %in% Type) {
    df <- cbind(Mu1_old, Mu2)
    colnames(df) <- c(GroupLabel1, GroupLabel2)
    mdf <- reshape2::melt(df)
    ggplot2::ggplot(mdf, aes(x = Var2, y = value)) +
      labs(title = "Before correction", y = "Mean expression", x = NULL) +
      geom_boxplot() +
      scale_y_log10()

    df <- cbind(Mu1, Mu2)
    colnames(df) <- c(GroupLabel1, GroupLabel2)
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

    g1 <- ggplot(mapping = aes(x = log2(MuBase_old), y = MedianTau_old)) +
      # geom_hex() + 
      ggplot2::geom_point() +
      labs(
        x = "Mean expresssion (log2)",
        y = paste(
          "Log2 fold change", GroupLabel1,
          "vs", GroupLabel2),
        title = "Before correction") +
      viridis::scale_fill_viridis() +
      ggplot2::geom_hline(yintercept = 0, lty = 2) +
      ggplot2::geom_hline(yintercept = log2(OffsetEst), lty = 1, col = "red")

    g2 <- ggplot2::ggplot(mapping = ggplot2::aes(x = log2(MuBase), y = MedianTau)) +
      # ggplot2::geom_hex() + 
      ggplot2::geom_point() +
      ggplot2::labs(
        x = "Mean expresssion (log2)",
        y = paste(
          "Log2 fold change", GroupLabel1,
          "vs", GroupLabel2),
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











BASiCS_plotDE <- function() {
  if (Search) {
    message("Plots to follow: \n",
            "1. EFDR/EFNR control plots \n",
            "2. MA plots \n",
            "3. Volcano plots \n")
  }
  else {
    message("Plots to follow: \n",
            "1. MA plots \n",
            "2. Volcano plots \n")
  }

  par(ask = TRUE)

  if (Search) {
    if (!is.null(Chain1@parameters$epsilon)) {
      par(mfrow = c(1, 3))
    }
    else {
      par(mfrow = c(1, 2))
    }
    ProbThresholds <- seq(0.5, 0.9995, by = 0.00025)
    plot(ProbThresholds, AuxMean$EFDRgrid,
         type = "l", lty = 1, bty = "n",
         ylab = "Error rate", xlab = "Probability threshold",
         ylim = c(0, 1), main = "Differential mean")
    lines(ProbThresholds, AuxMean$EFNRgrid, lty = 2)
    abline(h = EFDR_M, col = "blue", lwd = 2, lty = 1)
    abline(v = OptThresholdM[1], col = "red", lwd = 2, lty = 1)
    legend("top", c("EFDR", "EFNR", "Target EFDR"), lty = c(1, 2, 1),
           col = c("black", "black", "blue"), bty = "n")
    plot(ProbThresholds, AuxDisp$EFDRgrid,
         type = "l", lty = 1, bty = "n",
         ylab = "Error rate", xlab = "Probability threshold",
         ylim = c(0, 1), main = "Differential dispersion")
    lines(ProbThresholds, AuxDisp$EFNRgrid, lty = 2)
    abline(h = EFDR_D, col = "blue", lwd = 2, lty = 1)
    abline(v = OptThresholdD[1], col = "red", lwd = 2, lty = 1)
    legend("top", c("EFDR", "EFNR", "Target EFDR"), lty = c(1, 2, 1),
           col = c("black", "black", "blue"), bty = "n")
    if(!is.null(Chain1@parameters$epsilon)){
      plot(ProbThresholds, AuxResDisp$EFDRgrid, type = "l", lty = 1, bty = "n",
           ylab = "Error rate", xlab = "Probability threshold",
           ylim = c(0, 1), main = "Differential residual dispersion")
      lines(ProbThresholds, AuxResDisp$EFNRgrid, lty = 2)
      abline(h = EFDR_R, col = "blue", lwd = 2, lty = 1)
      abline(v = OptThresholdE[1], col = "red", lwd = 2, lty = 1)
      legend("top", c("EFDR", "EFNR", "Target EFDR"), lty = c(1, 2, 1),
             col = c("black", "black", "blue"), bty = "n")
    }
  }

  # MA plots
  if(!is.null(Chain1@parameters$epsilon)){
    par(mfrow = c(1, 3))
  }
  else {
    par(mfrow = c(1, 2))
  }
  with(TableMean,
       graphics::smoothScatter(log2(MeanOverall), MeanLog2FC,
                               bty = "n",
                               xlab = "Mean expresssion (log2)",
                               ylab = paste("Log2 fold change",
                                            GroupLabel1, "vs",
                                            GroupLabel2),
                               main = "Differential mean"))
  with(TableMean[!(TableMean$ResultDiffMean %in%
                   c("ExcludedByUser", "NoDiff")), ],
      points(log2(MeanOverall), MeanLog2FC, pch = 16, col = "red"))
  abline(h = c(-EpsilonM, EpsilonM), lty = 2)
  with(TableDisp,
       graphics::smoothScatter(log2(MeanOverall), DispLog2FC,
                               bty = "n",
                               xlab = "Mean expresssion (log2)",
                               ylab = paste("Log2 fold change",
                                            GroupLabel1, "vs",
                                            GroupLabel2),
                               main = "Differential dispersion"))
  with(TableDisp[!(TableDisp$ResultDiffDisp %in%
                     c("ExcludedFromTesting", "ExcludedByUser", "NoDiff")), ],
      points(log2(MeanOverall), DispLog2FC, pch = 16, col = "red"))
  abline(h = c(-EpsilonD, EpsilonD), lty = 2)

  if (!is.null(Chain1@parameters$epsilon)) {

    with(TableResDisp[TableResDisp$ResultDiffResDisp != "ExcludedFromTesting", ],
         graphics::smoothScatter(log2(MeanOverall), ResDispDistance,
                                 bty = "n",
                                 xlab = "Mean expresssion (log2)",
                                 ylab = paste("Difference",
                                              GroupLabel1, "vs",
                                              GroupLabel2),
                                 main = "Differential residual dispersion"))
    with(TableResDisp[!(TableResDisp$ResultDiffResDisp %in%
                       c("ExcludedFromTesting", "ExcludedByUser", "NoDiff")), ],
         points(log2(MeanOverall), ResDispDistance, pch = 16, col = "red"))
    abline(h = c(-EpsilonR, EpsilonR), lty = 2)
  }

  # Volcano plots
  if (!is.null(Chain1@parameters$epsilon)) {
    par(mfrow = c(1, 3))
  }
  else {
    par(mfrow = c(1, 2))
  }
  with(TableMean,
       graphics::smoothScatter(MeanLog2FC, ProbDiffMean,
                               bty = "n", ylim = c(0, 1),
                               ylab = "Posterior probability",
                               xlab = paste("Log2 fold change",
                                            GroupLabel1, "vs",
                                            GroupLabel2),
                               main = "Differential mean test"))
  with(TableMean[!(TableMean$ResultDiffMean %in%
                   c("ExcludedByUser", "NoDiff")), ],
       points(MeanLog2FC, ProbDiffMean, pch = 16, col = "red"))
  abline(v = c(-EpsilonM, EpsilonM), lty = 2)
  with(TableDisp,
       graphics::smoothScatter(DispLog2FC, ProbDiffDisp,
                               bty = "n", ylim = c(0, 1),
                               ylab = "Posterior probability",
                               xlab = paste("Log2 fold change",
                                            GroupLabel1, "vs",
                                            GroupLabel2),
                               main = "Differential dispersion test"))
  with(TableDisp[!(TableDisp$ResultDiffDisp %in%
                     c("ExcludedFromTesting", "ExcludedByUser", "NoDiff")), ],
       points(DispLog2FC, ProbDiffDisp, pch = 16, col = "red"))
  abline(v = c(-EpsilonD, EpsilonD), lty = 2)

  if (!is.null(Chain1@parameters$epsilon)) {

    with(TableResDisp[TableResDisp$ResultDiffResDisp != "ExcludedFromTesting", ],
         graphics::smoothScatter(ResDispDistance, ProbDiffResDisp,
                                 bty = "n", ylim = c(0, 1),
                                 ylab = "Posterior probability",
                                 xlab = paste("Difference",
                                              GroupLabel1, "vs",
                                              GroupLabel2),
                                 main = "Differential residual dispersion test"))
    with(TableResDisp[!(TableResDisp$ResultDiffResDisp %in%
                       c("ExcludedFromTesting", "ExcludedByUser", "NoDiff")), ],
         points(ResDispDistance, ProbDiffResDisp, pch = 16, col = "red"))
    abline(v = c(-EpsilonR, EpsilonR), lty = 2)
  }

  par(ask = FALSE)
}
