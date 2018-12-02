run_BASiCS_MCMC <- function(...) {
  capture.output(out <- BASiCS_MCMC(...))
  out
}
