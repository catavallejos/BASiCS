# BASiCS
Bayesian Analysis of Single-Cell Sequencing Data is an integrated Bayesian hierarchical model where:

- Cell-specific normalization constants are estimated as part of the model parameters.
- Technical variability is quantified based on spike-in genes that are artificially introduced to each analysed cells lysate.
- The total variability of the expression counts is decomposed into technical and biological components.

## Installation

BASiCS will be submitted to BioConductor. In the meantime, the development version can be installed from GitHub:

```R
# install.packages("devtools")
devtools::install_github("catavallejos/BASiCS", build_vignettes = TRUE)
```

This installation might fail if some of the dependency libraries are not yet installed. If so, please run the following lines and repeat the installation. 

```{r dependencies}
#library(devtools)
#source("http://bioconductor.org/biocLite.R")
#biocLite("BiocGenerics")
#biocLite("scran")
#install.packages("Rcpp")
```

After a successful installation, the BASiCS library can be loaded using[^footnoteInstall] 

```{r load_packages}
library(BASiCS)
```

[^footnoteInstall]: The warning `"No methods found in "BiocGenerics""` might appear. Please ignore. `BASiCS` imports some of the generic functions provided by `BiocGenerics` that do not have any methods attached.

### Troubleshooting

#### `gfortran` error

If you are unable to install BASiCS in a Mac OS X system, <a href="http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/"> you might require an additional gfortran library. </a> 

- If you are running R versions 3.0.0 - R 3.3.0, this can be installed in the Mac terminal using the following command:

```{}
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```
- If you are running R 3.4.0, you should install the official gfortran 6.1.0 binary for OS X El Capitan (10.11) from https://gcc.gnu.org/wiki/GFortranBinaries. 

#### `pandoc` error

If you are not using RStudio, you might need to manually install `pandoc` in order to create the vignette provided for BASiCS

Additional installation issues/solutions can be found here: https://github.com/catavallejos/BASiCS/issues/9

## How to use BASiCS?

- To detect highly and lowly variable genes within a populations of cells: please refer to the vignette

```R
vignette('BASiCSIntro')
```

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

Catalina A. Vallejos (cvallejos 'at' turing.ac.uk)

## Collaborators

- <a href="https://github.com/horta"> Danilo Horta </a>
- <a href="https://github.com/nilseling"> Nils Eling </a>

## Acknowledgements

We thank several members of the Marioni laboratory (EMBL-EBI; CRUK-CI) for support and discussions throughout the development of this R library. In particular, we are grateful to Aaron Lun (@LTLA, CRUK-CI) for advise and support during the preparation the Bioconductor submission. We also acknowledge feedback and contributions from (Github aliases provided within parenthesis): Ben Dulken (@bdulken), Chang Xu (@xuchang116), Danilo Horta (@Horta), Dmitriy Zhukov (@dvzhukov), Jens Preußner (@jenzopr), Joanna Dreux (@Joannacodes), Kevin Rue-Albrecht (@kevinrue), Luke Zappia (@lazappi) Simon Anders (@s-andrews), Yongchao Ge and Yuan Cao (@yuancao90), among others. 

This work has been supported by the MRC Biostatistics Unit (MRC grant no. MRC_MC_UP_0801/1; Catalina Vallejos and Sylvia Richardson), EMBL European Bioinformatics Institute (core European Molecular Biology Laboratory funding; Catalina Vallejos, Nils Eling and John Marioni), CRUK Cambridge Institute (core CRUK funding; John Marioni) and The Alan Turing Institute (EPSRC grant no. EP/N510129/1; Catalina Vallejos). 

## References

- <a href="http://dx.doi.org/10.1371/journal.pcbi.1004333">Vallejos <em>et al.</em> (2015). BASiCS: Bayesian Analysis of Single-Cell Sequencing Data </a>
- <a href="http://dx.doi.org/10.1186/s13059-016-0930-3">Vallejos <em>et al.</em> (2016). Beyond comparisons of means: understanding changes in gene expression at the single cell level</a>
