\name{ParamMouse}
\alias{ParamMouse}
\docType{data}
\title{
Estimated parameters for Mouse ESC dataset. 
}
\description{
Posterior medians of BASiCS model parameters for a subset of the Mouse ESC dataset provided in Islam et al (2014), as analysed in  Vallejos et al (2015). This dataset also includes IDs of the genes and cells used for the analysis.

For technical information about the experiment, please refer to Islam et al (2014). 
}
\usage{data(ParamMouse)}
\format{
  The format is:
  
List of 8 elements 

 \code{Mu}   : Gene-specific expression levels, defined as true input molecules in case of technical genes 
 (vector of length 7941, including 7895 intrinsinc genes and 46 spike-in genes).
 
 \code{Delta}   : Gene-specific biological cell-to-cell variability hyper-parameters, biological genes only 
 (vector of length 7895). 
 
 \code{Phi}   : Cell-specific mRNA content normalising constants 
 (vector of length 41).
 
 \code{S}   : Cell-specific capture efficiency (or amplification biases if not using UMI based counts) normalising constants
 (vector of length 41). 
 
 \code{Nu}   : Cell-specific random effects related to the technical noise component
 (vector of length 41).
 
 \code{Theta}   : Technical variability hyper-parameter
 (a single positive number). 
 
 \code{CellsIDs}   : Identifyiers of the analysed cells, i.e. those cells who passed the quality control criterion adopted by Islam et al (vector of names with length 41).
 
 \code{GeneIDs}   : Identifyiers of the analysed genes, i.e. those genes that have (on average) at least 1 count across the 41 analysed cells (vector of names with length 7941, including 7895 intrinsinc genes and 46 spike-in genes).
 
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
\href{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46980}{Series GSE46980 (File: GSE46980_CombinedMoleculeCounts.tab.gz)}
}
\references{
Islam, Zeisel, Joost, La Manno, Zajac, Kasper, Lonnerberg and Linnarson (2014). Quantitative single-cell RNA-seq with unique molecular identifiers. Nature Methods, Vol 11, No 2 (\href{http://www.nature.com/nmeth/journal/v11/n2/full/nmeth.2772.html}{doi:10.1038/nmeth.2772}).

Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. 
}
\examples{
data(ParamMouse)
}
\keyword{datasets}
