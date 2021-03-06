% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/newBASiCS_Data.R
\name{newBASiCS_Data}
\alias{newBASiCS_Data}
\title{Creates a SingleCellExperiment object from a matrix of expression
counts and experimental information about spike-in genes}
\usage{
newBASiCS_Data(
  Counts,
  Tech = rep(FALSE, nrow(Counts)),
  SpikeInfo = NULL,
  BatchInfo = NULL,
  SpikeType = "ERCC"
)
}
\arguments{
\item{Counts}{Matrix of dimensions \code{q} times \code{n} whose elements
contain the expression counts to be analyses
(including biological and technical spike-in genes). Gene names must be
stored as \code{rownames(Counts)}.}

\item{Tech}{Logical vector of length \code{q}. If \code{Tech = FALSE} the
gene is biological; otherwise the gene is spike-in. Defaul value:
\code{Tech = rep(FALSE, nrow(Counts))}.}

\item{SpikeInfo}{\code{data.frame} whose first and second columns contain
the gene names assigned to the spike-in genes (they must match the ones in
\code{rownames(Counts)}) and the associated input number of molecules,
respectively. If \code{SpikeInfo = NULL}, only the horizontal integration
implementation (no spikes) can be run. Default value:
\code{SpikeInfo = NULL}.}

\item{BatchInfo}{Vector of length \code{n} whose elements indicate batch
information. Not required if a single batch is present on the data.
Default value: \code{BatchInfo = NULL}.}

\item{SpikeType}{Character to indicate what type of spike-ins are in use.
Default value: \code{SpikeType = "ERCC"} (parameter is no longer used).}
}
\value{
An object of class \code{\linkS4class{SingleCellExperiment}}.
}
\description{
\code{newBASiCS_Data} creates a
\code{\linkS4class{SingleCellExperiment}} object from a matrix of expression
counts and experimental information about spike-in genes.
}
\seealso{
\code{\linkS4class{SingleCellExperiment}}
}
\author{
Catalina A. Vallejos \email{cnvallej@uc.cl}

Nils Eling \email{eling@ebi.ac.uk}
}
