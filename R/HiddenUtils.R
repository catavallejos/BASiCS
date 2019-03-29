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
                                         "geweke.diag"), 
                             na.rm = FALSE) {

  Measure <- match.arg(Measure)
  MeasureFun <- match.fun(Measure)
  mat <- HiddenGetParam(object, Param)
  if (na.rm) {
    mat <- mat[, !apply(mat, 2, function(col) any(is.na(col)))]
  }
  metric <- MeasureFun(coda::mcmc(mat))
  if (Measure == "geweke.diag") {
    metric <- metric$z
  }
  metric
}
