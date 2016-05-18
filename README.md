# BASiCS
Bayesian Analysis of Single-Cell Sequencing Data is an integrated Bayesian hierarchical model where:

- Cell-specific normalization constants are estimated as part of the model parameters.
- Technical variability is quantified based on spike-in genes that are artificially introduced to each analysed cells lysate.
- The total variability of the expression counts is decomposed into technical and biological components.

## Installation

BASiCS will be submitted to BioConductor. In the meantime, the development version can be install from GitHub:

```R
# install.packages("devtools")
devtools::install_github("catavallejos/BASiCS")
```

This installation might fail if some of the dependency libraries are not yet installed. If so, please run the following lines and repeat the installation. 

```{r dependencies}
#library(devtools)
#source("http://bioconductor.org/biocLite.R")
#biocLite("BiocGenerics")
#install.packages("Rcpp")
```

After a successful installation, the BASiCS library can be loaded using[^footnoteInstall] 

```{r load_packages}
library(BASiCS)
```

[^footnoteInstall]: The warning `"No methods found in "BiocGenerics""` might appear. Please ignore. `BASiCS` imports some of the generic functions provided by `BiocGenerics` that do not have any methods attached.

## How to use BASiCS?

- To detect highly and lowly variable genes within a populations of cells: please refer to the <a href="https://github.com/catavallejos/BASiCS/blob/master/vignettes/BASiCSIntro.Rmd">BASiCS vignette</a>

- To detect changes whose expression changes between 2 or more populations of cells (mean and over-dispersion), please refer to the supplementary material related to <a href="http://dx.doi.org/10.1186/s13059-016-0930-3">Vallejos <em>et al.</em>, 2016</a>

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

## Collaborators

- <a href="https://github.com/horta"> Danilo Horta </a>

## References

- <a href="http://dx.doi.org/10.1371/journal.pcbi.1004333">Vallejos <em>et al.</em> (2015). BASiCS: Bayesian Analysis of Single-Cell Sequencing Data </a>
- <a href="http://dx.doi.org/10.1186/s13059-016-0930-3">Vallejos <em>et al.</em> (2016). Beyond comparisons of means: understanding changes in gene expression at the single cell level</a>
