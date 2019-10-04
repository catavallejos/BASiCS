HiddenScaleName <- function(Measure = c("effectiveSize",
                                        "geweke.diag"),
                            Param = NULL) {
  Measure <- match.arg(Measure)
  measure_name <- switch(Measure, 
    effectiveSize = "Effective sample size",
    geweke.diag = "Geweke diagnostic"
  )
  if (!is.null(Param)) {
    measure_name <- paste0(measure_name, ":\n", Param)
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
    if (!ncol(mat)) {
      stop(paste("No non-NA samples for", Param))
    }
  }
  metric <- MeasureFun(coda::mcmc(mat))
  if (Measure == "geweke.diag") {
    metric <- metric$z
  }
  metric
}

HiddenGetParam <- function(object, Param = "mu") {
  if (is.null(Param) || 
      is.na(Param) || 
      length(Param) > 1 ||
      !(Param %in% names(object@parameters))) {
    stop("'Param' argument is invalid")
  }
  object@parameters[[Param]]
}

HiddenCheckValidCombination <- function(...) {
  Params <- list(...)
  Check1 <- vapply(Params, 
                   FUN = function(x) is.null(x) || x %in% HiddenGeneParams(),
                   FUN.VALUE = TRUE)
  Check2 <- vapply(Params, 
                   FUN = function(x) is.null(x) || x  %in% HiddenCellParams(),
                   FUN.VALUE = TRUE)
  
  if (!(all(Check1) || all(Check2))) {
    stop(paste("Invalid combination of parameters:",
               paste(list(...), collapse = ", "), " \n"))
  } 
}

HiddenGeneParams <- function() c("mu", "delta", "epsilon")
HiddenCellParams <- function() c("s", "phi", "nu")

NClassFD2D <- function(x, y) {
  max(nclass.FD(x), nclass.FD(y))
}

MeasureName <- function(measure) {
  switch(measure, 
    "Mean" = "mean expression", 
    "Disp" = "over dispersion",
    "ResDisp" = "residual over dispersion")
}

DistanceName <- function(measure) {
  switch(measure, 
    "ResDisp" = "distance",
    "fold change")
}

LogDistanceName <- function(measure) {
  switch(measure, 
    "ResDisp" = "distance",
    "log2(fold change)")
}

DistanceVar <- function(measure) {
  switch(measure, 
    "ResDisp" = "Distance",
    "Log2FC")
}


cap <- function(s) {
  sub("([[:alpha:]])([[:alpha:]]+)", "\\U\\1\\L\\2", s, perl = TRUE)
}
NSamples <- function(Chain) {
  nrow(Chain@parameters[[1]])
}
