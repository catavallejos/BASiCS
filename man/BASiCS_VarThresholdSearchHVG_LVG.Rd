% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BASiCS_VarThresholdSearchHVG_LVG.R
\name{BASiCS_VarThresholdSearchHVG}
\alias{BASiCS_VarThresholdSearchHVG}
\alias{BASiCS_VarThresholdSearchVG}
\alias{BASiCS_VarThresholdSearchHVG_LVG}
\alias{BASiCS_VarThresholdSearchLVG}
\title{Detection method for highly and lowly variable genes using
a grid of variance contribution thresholds}
\usage{
BASiCS_VarThresholdSearchVG(
  Chain,
  Task = c("HVG", "LVG"),
  VarThresholdsGrid,
  EFDR = 0.1,
  Progress = TRUE
)

BASiCS_VarThresholdSearchHVG(...)

BASiCS_VarThresholdSearchLVG(...)
}
\arguments{
\item{Chain}{an object of class \code{\linkS4class{BASiCS_Chain}}}

\item{Task}{See \code{?BASiCS_DetectVG}.}

\item{VarThresholdsGrid}{Grid of values for the variance contribution
threshold (they must be contained in (0,1))}

\item{EFDR}{Target for expected false discovery rate related to
HVG/LVG detection. Default: \code{EFDR = 0.10}.}

\item{Progress}{If \code{Progress = TRUE}, partial output is
printed in the console. Default: \code{Progress = TRUE}.}

\item{...}{Passed to methods.}
}
\value{
\describe{
\item{\code{BASiCS_VarThresholdSearchHVG}}{A table displaying the results of
      highly variable genes detection for different variance
      contribution thresholds.}
\item{\code{BASiCS_VarThresholdSearchLVG}}{A table displaying the results of
      lowly variable genes detection for different variance
      contribution thresholds.}
}
}
\description{
Detection method for highly and lowly variable genes
using a grid of variance contribution thresholds. Only used when 
HVG/LVG are found based on the variance decomposition.
}
\details{
See vignette
}
\examples{

data(ChainSC)

BASiCS_VarThresholdSearchHVG(ChainSC,
                             VarThresholdsGrid = seq(0.55,0.65,by=0.01),
                             EFDR = 0.10)
BASiCS_VarThresholdSearchLVG(ChainSC,
                             VarThresholdsGrid = seq(0.35,0.45,by=0.01),
                             EFDR = 0.10)

}
\references{
Vallejos, Marioni and Richardson (2015). PLoS Computational Biology.
}
\seealso{
\code{\linkS4class{BASiCS_Chain}}
}
\author{
Catalina A. Vallejos \email{cnvallej@uc.cl}
}
