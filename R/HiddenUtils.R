## Slightly optimised version of 
## coda::effectiveSize 
## https://cran.r-project.org/web/packages/coda/
#' @importFrom matrixStats colVars
#' @importFrom stats ar setNames
ess <- function(x) {
  vars <- matrixStats::colVars(x)
  spec <- numeric(ncol(x))
  has_var <- vars != 0
  if (any(has_var, na.rm = TRUE)) {
    spec[which(has_var)] <- apply(x[, which(has_var)],
      2,
      function(y) {
        a <- ar(y, aic = TRUE)
        a$var.pred / (1 - sum(a$ar)) ^ 2
      }
    )    
  }
  setNames(ifelse(spec == 0, 0, nrow(x) * vars / spec), colnames(x))
}

HiddenScaleName <- function(Measure = c("ess",
                                        "geweke.diag"),
                            Param = NULL) {
  Measure <- match.arg(Measure)
  measure_name <- switch(Measure, 
    ess = "Effective sample size",
    geweke.diag = "Geweke diagnostic"
  )
  if (!is.null(Param)) {
    measure_name <- paste0(measure_name, ": ", Param)
  }
  measure_name
}

HiddenGetMeasure <- function(object, 
                             Param,
                             Measure = c("ess",
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
    "ResDisp" = "residual over dispersion"
  )
}

DistanceName <- function(measure) {
  switch(measure,
    "ResDisp" = "distance",
    "fold change")
}

LogDistanceName <- function(measure) {
  switch(measure,
    "ResDisp" = "distance",
    "log2(fold change)"
  )
}

DistanceVar <- function(measure) {
  switch(measure,
    "ResDisp" = "Distance",
    "Log2FC"
  )
}

cap <- function(s) {
  sub("([[:alpha:]])([[:alpha:]]+)", "\\U\\1\\L\\2", s, perl = TRUE)
}

NSamples <- function(Chain) {
  nrow(Chain@parameters[[1]])
}
