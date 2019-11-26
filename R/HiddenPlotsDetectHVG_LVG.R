HiddenPlot1DetectHVG_LVG <- function(ProbThresholds, EFDRgrid, EFNRgrid, EFDR) {
  ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(ProbThresholds, EFDRgrid, color = "EFDR")) +
    ggplot2::geom_line(ggplot2::aes(ProbThresholds, EFNRgrid, color = "EFNR")) +
    ggplot2::geom_hline(ggplot2::aes(color = "Target EFDR", yintercept = EFDR)) +
    ggplot2::labs(x = "Probability threshold", y = "Error rate") +
    ggplot2::ylim(c(0, 1))
}


HiddenPlot2DetectHVG_LVG <- function(Task,
                                     Mu,
                                     Prob,
                                     OptThreshold,
                                     Hits,
                                     ylim = c(0, 1),
                                     xlim = c(min(Mu), max(Mu)),
                                     cex = 1.5,
                                     pch = 16,
                                     col = 8,
                                     bty = "n",
                                     cex.lab = 1,
                                     cex.axis = 1,
                                     cex.main = 1,
                                     xlab = "Mean expression",
                                     ylab = paste(Task, "probability"),
                                     main = "") {

  col <- grDevices::rgb(grDevices::col2rgb(col)[1],
                        grDevices::col2rgb(col)[2],
                        grDevices::col2rgb(col)[3], 50,
                        maxColorValue = 255)

  df <- data.frame(Mu, Prob)
  ggplot2::ggplot(df, ggplot2::aes_string(x = "Mu", y = "Prob")) +
    ggplot2::geom_point(
      aes(color = ifelse(Hits, Task, "Other")), 
      pch = pch, cex = cex) +
    ggplot2::scale_color_brewer(palette = "Set1", name = "") +
    ggplot2::geom_hline(yintercept = OptThreshold[[1]], lty = 2, col = "black") +
    # ggplot2::geom_point(data = df[Hits, ], pch = pch, col = "red", cex = cex) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      title = main
    )
}
