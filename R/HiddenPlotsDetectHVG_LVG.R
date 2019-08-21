HiddenPlot1DetectHVG_LVG <- function(ProbThresholds, EFDRgrid, EFNRgrid, EFDR) {
  plot(ProbThresholds, EFDRgrid,
       type = "l", lty = 1, bty = "n", lwd = 2,
       ylab = "Error rate", xlab = "Probability threshold",
       ylim = c(0, 1))
  lines(ProbThresholds, EFNRgrid, lty = 2, lwd = 2)
  abline(h = EFDR, col = "blue", lwd = 2, lty = 1)
  legend("topleft", c("EFDR", "EFNR", "Target EFDR"),
         lty = c(1, 2, 1), col = c("black", "black", "blue"),
         bty = "n", lwd = 2)
}
# HiddenPlot2DetectHVG_LVG <- function(Task,
#                                      Mu,
#                                      Prob,
#                                      OptThreshold,
#                                      Hits,
#                                      ylim = c(0, 1),
#                                      xlim = c(min(Mu), max(Mu)),
#                                      cex = 1.5,
#                                      pch = 16,
#                                      col = 8,
#                                      bty = "n",
#                                      cex.lab = 1,
#                                      cex.axis = 1,
#                                      cex.main = 1,
#                                      xlab = "Mean expression",
#                                      ylab = paste(Task, "probability"),
#                                      main = "") {

#   col <- grDevices::rgb(grDevices::col2rgb(col)[1],
#                         grDevices::col2rgb(col)[2],
#                         grDevices::col2rgb(col)[3], 50,
#                         maxColorValue = 255)

#   plot(Mu, Prob, log = "x", pch = pch,
#        ylim = ylim, xlim = xlim, col = col,
#        cex = cex, bty = bty, cex.lab = cex.lab,
#        cex.axis = cex.axis, cex.main = cex.main,
#        xlab = xlab, ylab = ylab, main = main)
#   abline(h = OptThreshold[1], lty = 2, col = "black")
#   points(Mu[Hits], Prob[Hits], pch = pch,
#          col = grDevices::rgb(grDevices::col2rgb("red")[1],
#                               grDevices::col2rgb("red")[2],
#                               grDevices::col2rgb("red")[3],
#                               50, maxColorValue = 255), cex = cex)
# }



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
  ggplot(df, aes_string(x = "Mu", y = "Prob")) +
    geom_point(pch = pch, cex = cex) +
    geom_hline(yintercept = OptThreshold[[1]], lty = 2, col = "black") +
    geom_point(data = df[Hits, ], pch = pch, col = "red", cex = cex) +
    labs(
      x = xlab,
      y = ylab,
      title = main
    )
}
