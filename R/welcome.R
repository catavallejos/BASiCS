# ==========================================================================
# package initialization
# ==========================================================================
.onAttach = function(libname, pkgname) {
  #msg = "BASiCS (v0.5.7 onwards) has been extended to cover differential expression analyses (mean and over-dispersion), as described in http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0930-3. New vignette file is in preparation. In the meantime, please refer to the supplementary files provided in the link for the code required to use the new functionalities. If you have any doubts, please contact cnvallej@uc.cl. "
  msg = "BASiCS now (v1.0.0 onwards) utilizes the SummarizedExperiment class to collect the dataset information: expression counts, information on spike-ins and batch information. For more information use vignette('BASiCS'). \n 
  In preparation to the Bioconductor submission, multiple functions have been modified. \n
  We aimed to minimise changes that affect usage syntax, but some were required to improve functionality. \n
  For a summary of the changes, please visit: https://github.com/catavallejos/BASiCS/wiki."
  msg = strwrap(msg, exdent=4, indent=4)
  packageStartupMessage(paste(msg, collapse="\n"), appendLF=TRUE)
}

# Add message describing which functions changed
# And message describing which functions/classes were removed
