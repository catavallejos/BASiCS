# ==============================================================
# package initialization
# ==============================================================
.onAttach = function(libname, pkgname) {
  msg = "BASiCS now (v1.0.0 onwards) utilizes the SingleCellExperiment class \n 
  to collect the dataset information: expression counts, information on spike-ins \n 
  and batch information. For more information use vignette('BASiCS'). \n 
  In preparation to the Bioconductor submission, multiple functions have been modified. \n
  We aimed to minimise changes that affect usage syntax, \n 
  but some were required to improve functionality. \n
  For a summary of the changes, please visit: https://github.com/catavallejos/BASiCS/wiki."
  msg = strwrap(msg, exdent=4, indent=4)
  packageStartupMessage(paste(msg, collapse="\n"), appendLF=TRUE)
}