HiddenScaleName <- function(Measure = c("effectiveSize",
                                  "geweke.diag"),
                      Param = NULL) {
  Measure <- match.arg(Measure)
  measure_name <- switch(Measure, 
    effectiveSize = "Effective sample size",
    geweke.diag = "Geweke diagnostic"
  )
  if (!is.null(Param)) {
    measure_name <- paste0(measure_name, ": ", Param)
  }
  measure_name
}

HiddenGetMeasure <- function(object, 
                       Param,
                       Measure = c("effectiveSize",
                                   "geweke.diag")) {
  Measure <- match.arg(Measure)
  MeasureFun <- match.fun(Measure)
  metric <- MeasureFun(coda::mcmc(HiddenGetParam(object, Param)))
  if (Measure == "geweke.diag") {
    metric <- metric$z
  }
  metric
}
