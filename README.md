# BASiCS
Bayesian Analysis of Single-Cell Sequencing Data is an integrated Bayesian hierarchical model where:

- Cell-specific normalization constants are estimated as part of the model parameters.
- Technical variability is quantified based on spike-in genes that are artificially introduced to each analysed cells lysate.
- The total variability of the expression counts is decomposed into technical and biological components.

## Installation

Get the released version from CRAN:

```R
install.packages("BASiCS")
```

Or install the development version from GitHub:

```R
# install.packages("devtools")
# source("https://bioconductor.org/biocLite.R")
# biocLite("BiocGenerics")
devtools::install_github("catavallejos/BASiCS")
```

## Quick start

TODO: a quick start for BASiCS. Like vignette("some-stuff")

## Development

### Sanity check

In the folder of the development version from GitHub:

```
R CMD check .
```

### Build from source

In the parent folder of the development version from GitHub:

```
R CMD build BASiCS
```

And run then

```R
install.packages('BASiCS_x.x.x.tar.gz', repos=NULL)
```

on the generated file.

## Author

Catalina A. Vallejos <catalina.vallejos@mrc-bsu.cam.ac.uk>

## References

- [BASiCS: Bayesian Analysis of Single-Cell Sequencing Data](http://www.ncbi.nlm.nih.gov/pubmed/26107944)
