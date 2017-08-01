# BASiCS

Bayesian Analysis of Single-Cell Sequencing Data is an integrated Bayesian hierarchical model to perform statistical analyses of single-cell RNA sequencing datasets in the context of **supervised** experiments (where the groups of cells of interest are known a priori, e.g. experimental conditions or cell types). 

In BASiCS:

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

```R
#library(devtools)
#source("http://bioconductor.org/biocLite.R")
#biocLite("BiocGenerics")
#biocLite("scran")
#install.packages("Rcpp")
```

After a successful installation, the BASiCS library can be loaded using[^footnoteInstall] 

```R
library(BASiCS)
```

[^footnoteInstall]: The warning `"No methods found in "BiocGenerics""` might appear. Please ignore. `BASiCS` imports some of the generic functions provided by `BiocGenerics` that do not have any methods attached.

## Installation troubleshooting

A summary of the installation errors that have been reported for BASiCS is provided [here](https://github.com/catavallejos/BASiCS/wiki/7.-Installation-troubleshooting). If you encounter any additional issues, **please let us know so that we can update this information**.

## How to use BASiCS?

BASiCS includes a vignette where its usage is illutrated. To access the vignette, please use:

```R
vignette('BASiCS')
```

Individual topics are summarized in the BASiCS wiki:

- [Quick start](https://github.com/catavallejos/BASiCS/wiki/1.-Quick-start)

- [Input preparation](https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation)

- [MCMC convergence](https://github.com/catavallejos/BASiCS/wiki/3.-MCMC-convergence)

- [Posterior summary](https://github.com/catavallejos/BASiCS/wiki/4.-Posterior-summary)

- [HVL & LVG detection](https://github.com/catavallejos/BASiCS/wiki/5.-HVG-&-LVG-detection) for a single group of cells

- [Differential analysis](https://github.com/catavallejos/BASiCS/wiki/6.-Differential-analysis) between 2 groups of cells (mean and over-dispersion)


<!---- To detect changes whose expression changes between 2 or more populations of cells (mean and over-dispersion), please refer to the supplementary material related to <a href="http://dx.doi.org/10.1186/s13059-016-0930-3">Vallejos <em>et al.</em>, 2016</a> TODO: a quick start for BASiCS. Like vignette("some-stuff") ---> 

## Authors

- [Catalina Vallejos](https://sites.google.com/view/catalinavallejos) (cvallejos 'at' turing.ac.uk)
- [Nils Eling](https://github.com/nilseling)
- John Marioni
- Sylvia Richardson

## Acknowledgements

We thank several members of the Marioni laboratory (EMBL-EBI; CRUK-CI) for support and discussions throughout the development of this R library. In particular, we are grateful to Aaron Lun (@LTLA, CRUK-CI) for advise and support during the preparation the Bioconductor submission. 

We also acknowledge feedback and contributions from (Github aliases provided within parenthesis): Ben Dulken (@bdulken), Chang Xu (@xuchang116), Danilo Horta (@Horta), Dmitriy Zhukov (@dvzhukov), Jens Preußner (@jenzopr), Joanna Dreux (@Joannacodes), Kevin Rue-Albrecht (@kevinrue), Luke Zappia (@lazappi), Simon Anders (@s-andrews), Yongchao Ge and Yuan Cao (@yuancao90), among others. 

This work has been funded by the MRC Biostatistics Unit (MRC grant no. MRC_MC_UP_0801/1; Catalina Vallejos and Sylvia Richardson), EMBL European Bioinformatics Institute (core European Molecular Biology Laboratory funding; Catalina Vallejos, Nils Eling and John Marioni), CRUK Cambridge Institute (core CRUK funding; John Marioni) and The Alan Turing Institute (EPSRC grant no. EP/N510129/1; Catalina Vallejos). 

## References

- <a href="http://dx.doi.org/10.1371/journal.pcbi.1004333">Vallejos <em>et al.</em> (2015). PLoS Computational Biology. </a>
- <a href="http://dx.doi.org/10.1186/s13059-016-0930-3">Vallejos <em>et al.</em> (2016). Genome Biology. </a>
