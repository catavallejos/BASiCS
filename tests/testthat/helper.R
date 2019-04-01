run_MCMC <- function(...) {
    capture.output(chain <- BASiCS_MCMC(...))
    chain
}
