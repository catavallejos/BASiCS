[![Build Status](https://travis-ci.org/catavallejos/BASiCS.svg?branch=master)](https://travis-ci.org/catavallejos/BASiCS) 
![codecov.io Code Coverage](https://codecov.io/gh/catavallejos/BASiCS/branch/master/graph/badge.svg)

# BASiCS

Single-cell mRNA sequencing can uncover novel cell-to-cell heterogeneity in gene 
expression levels within seemingly homogeneous populations of cells. However, 
these experiments are prone to high levels of technical noise, creating new 
challenges for identifying genes that show genuine heterogeneous expression 
within the group of cells under study. 

BASiCS (**B**ayesian **A**nalysis of **Si**ngle-**C**ell **S**equencing data) is 
an integrated Bayesian hierarchical model that propagates statistical 
uncertainty by simultaneously performing data normalisation (global scaling), 
technical noise quantification and two types of **supervised** downstream
analyses: 

- **For a single group of cells** [1]: BASiCS provides a criterion to identify 
highly (and lowly) variable genes within the group. 

- **For two (or more) groups of cells** [2]: BASiCS allows the identification 
of differentially expressed genes between the groups. As in traditional 
differential expression tools, BASiCS can uncover changes in mean expression 
between the groups. Besides this, BASiCS can also uncover changes in 
*over-dispersion* --- a measure for the excess cell-to-cell variation that is 
observed after accounting for technical noise. This feature has led, 
for example, to novel insights in the context of immune cells across aging [3]. 
More recently, the BASiCS model has been extended to address the confounding
between mean and variability that is typically observed in scRNA-seq datasets.
This is achieved by introducing a *residual over-dispersion* parameter that 
is not confounded by mean expression [4]. 

In both cases, a probabilistic output is provided, with posterior probability 
thresholds calibrated through the expected false discovery rate (EFDR) [5].

The original implementation of BASiCS relies on the use of **spike-in genes** 
--- that are artificially introduced to each cell's lysate --- to perform these 
analyses. However, our latest work has extended the BASiCS model to datasets
in which spike-ins are not available (multiple *batches* are required) [4].


**Important: BASiCS has been designed in the context of supervised experiments where the groups of cells (e.g. experimental conditions, cell types) under study are known a priori (e.g. case-control studies). Therefore, we DO NOT advise the use of BASiCS in unsupervised settings where the aim is to uncover sub-populations of cells through clustering.**

For technical details, references are provided at the bottom of this document. 

## Installation

BASiCS is available in [Bioconductor](https://bioconductor.org/packages/BASiCS).
To install the current release use:

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("BASiCS")
```

Repeat using the [devel](https://bioconductor.org/developers/how-to/useDevel/) 
version of Bioconductor for the latest development version. 

Alternatively, the experimental version of BASiCS (this might be unstable)
can be installed from GitHub:

```R
# install.packages("devtools")
devtools::install_github("catavallejos/BASiCS", build_vignettes = TRUE)
```

This installation might fail if some of the dependency libraries are not yet 
installed. If so, please run the following lines and repeat the installation. 

```R
#library(devtools)
#if (!requireNamespace("BiocManager", quietly=TRUE))
    #install.packages("BiocManager")
#BiocManager::install("BiocGenerics")
#BiocManager::install("scran")
#install.packages("Rcpp")
```

After a successful installation, the BASiCS library can be 
loaded using[^footnoteInstall] 

```R
library(BASiCS)
```

[^footnoteInstall]: The warning `"No methods found in "BiocGenerics""` might 
appear. Please ignore. `BASiCS` imports some of the generic functions provided 
by `BiocGenerics` that do not have any methods attached.

## Installation troubleshooting

A summary of the installation errors that have been reported for BASiCS is 
provided [here](https://github.com/catavallejos/BASiCS/wiki/7.-Installation-troubleshooting). 
If you encounter any additional issues, **please let us know so that we can update this information**.

## How to use BASiCS?

BASiCS includes a vignette where its usage is illutrated. 
To access the vignette, please use:

```R
vignette('BASiCS')
```

Individual topics are summarized in the BASiCS wiki:

- [Quick start](https://github.com/catavallejos/BASiCS/wiki/1.-Quick-start)

- [Input preparation](https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation)

- [Running the MCMC](https://github.com/catavallejos/BASiCS/wiki/3.-Running-the-MCMC)

- [MCMC convergence](https://github.com/catavallejos/BASiCS/wiki/4.-MCMC-convergence)

- [Posterior summary](https://github.com/catavallejos/BASiCS/wiki/5.-Posterior-summary)

- [HVL & LVG detection](https://github.com/catavallejos/BASiCS/wiki/6.-HVG-&-LVG-detection) for a single group of cells

- [Differential analysis](https://github.com/catavallejos/BASiCS/wiki/7.-Differential-analysis) between 2 groups of cells (mean and over-dispersion)


<!---- To detect changes whose expression changes between 2 or more populations of cells (mean and over-dispersion), please refer to the supplementary material related to <a href="http://dx.doi.org/10.1186/s13059-016-0930-3">Vallejos <em>et al.</em>, 2016</a> TODO: a quick start for BASiCS. Like vignette("some-stuff") ---> 

## Authors

- [Catalina Vallejos](https://sites.google.com/view/catalinavallejos) (cvallejos 'at' turing.ac.uk)
- [Nils Eling](https://github.com/nilseling)
- John Marioni
- Sylvia Richardson

## Acknowledgements

We thank several members of the Marioni laboratory (EMBL-EBI; CRUK-CI) for 
support and discussions throughout the development of this R library. 
In particular, we are grateful to Aaron Lun (@LTLA, CRUK-CI) for advise and 
support during the preparation the Bioconductor submission. 

We also acknowledge feedback and/or contributions from (Github aliases provided 
within parenthesis): Alan O'Callaghan (@Alanocallaghan), Ben Dulken (@bdulken), 
Chang Xu (@xuchang116), Danilo Horta (@Horta), Dmitriy Zhukov (@dvzhukov), 
Jens Preußner (@jenzopr), Joanna Dreux (@Joannacodes), Kevin Rue-Albrecht 
(@kevinrue), Luke Zappia (@lazappi), Mike Morgan (@MikeDMorgan), Muad Abd El Hay 
(@Cumol), Nitesh Turaga (@nturaga), Simon Anders (@s-andrews), Yongchao Ge and 
Yuan Cao (@yuancao90), among others. 

This work has been funded by the MRC Biostatistics Unit (MRC grant no. 
MRC_MC_UP_0801/1; Catalina Vallejos and Sylvia Richardson), 
EMBL European Bioinformatics Institute (core European Molecular Biology 
Laboratory funding; Catalina Vallejos, Nils Eling and John Marioni), 
CRUK Cambridge Institute (core CRUK funding; John Marioni) and The Alan Turing 
Institute (EPSRC grant no. EP/N510129/1; Catalina Vallejos). 

## References

- [1] <a href="http://dx.doi.org/10.1371/journal.pcbi.1004333">Vallejos <em>et al.</em> (2015). PLoS Computational Biology. </a>
- [2] <a href="http://dx.doi.org/10.1186/s13059-016-0930-3">Vallejos <em>et al.</em> (2016). Genome Biology. </a>
- [3] <a href="http://science.sciencemag.org/content/355/6332/1433">Martinez-Jimenes <em>et al.</em> (2017). Science. </a>
- [4] <a href="https://www.cell.com/cell-systems/fulltext/S2405-4712(18)30278-3">Eling <em>et al.</em> (2018). Cell Systems. </a>
- [5] <a href="https://www.ncbi.nlm.nih.gov/pubmed/15054023">Newton <em>et al.</em> (2004). Biostatistics. </a>
