HiddenPlot1DetectHVG_LVG <- function(ProbThresholds, EFDRgrid, EFNRgrid, EFDR) 
{
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

HiddenPlot2DetectHVG_LVG <- function(args, Task, Mu, Prob, OptThreshold, Hits) 
{
  if ("ylim" %in% names(args)) { ylim <- args$ylim } 
  else { ylim <- c(0, 1) }
  if ("xlim" %in% names(args)) { xlim <- args$xlim } 
  else { xlim <- c(min(Mu), max(Mu)) }
  
  cex <- ifelse("cex" %in% names(args), args$cex, 1.5)
  pch <- ifelse("pch" %in% names(args), args$pch, 16)
  col <- ifelse("col" %in% names(args), args$col, 8)
  bty <- ifelse("bty" %in% names(args), args$bty, "n")
  cex.lab <- ifelse("cex.lab" %in% names(args), args$cex.lab, 1)
  cex.axis <- ifelse("cex.axis" %in% names(args), args$cex.axis, 1)
  cex.main <- ifelse("cex.main" %in% names(args), args$cex.main, 1)
  
  xlab <- ifelse("xlab" %in% names(args), args$xlab, "Mean expression")
  if (Task == "HVG") {
    ylab <- ifelse("ylab" %in% names(args), args$ylab, "HVG probability")
  } else {
    ylab <- ifelse("ylab" %in% names(args), args$ylab, "LVG probability")
  }
  main <- ifelse("main" %in% names(args), args$main, "")
  
  col <- grDevices::rgb(grDevices::col2rgb(col)[1], 
                        grDevices::col2rgb(col)[2], 
                        grDevices::col2rgb(col)[3], 50, 
                        maxColorValue = 255)
  
  plot(Mu, Prob, log = "x", pch = pch, 
       ylim = ylim, xlim = xlim, col = col, 
       cex = cex, bty = bty, cex.lab = cex.lab, 
       cex.axis = cex.axis, cex.main = cex.main,
       xlab = xlab, ylab = ylab, main = main)
  abline(h = OptThreshold[1], lty = 2, col = "black")
  points(Mu[Hits], Prob[Hits], pch = pch, 
         col = grDevices::rgb(grDevices::col2rgb("red")[1], 
                              grDevices::col2rgb("red")[2],
                              grDevices::col2rgb("red")[3], 
                              50, maxColorValue = 255), cex = cex)
}