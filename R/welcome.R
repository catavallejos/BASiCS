# ==============================================================
# package initialization
# ==============================================================
.onAttach = function(libname, pkgname) {
  msg = "BASiCS now (v0.99.0 onwards) utilizes the SingleCellExperiment class \n 
  to collect the dataset information: expression counts, information on \n 
  spike-ins and batch information. For more information use \n 
  vignette('BASiCS'). In preparation to the Bioconductor submission, multiple \n 
  functions have been modified. We aimed to minimise changes that affect \n 
  usage syntax, but some were required to improve functionality. For a \n 
  summary of the changes, please visit: \n
  https://github.com/catavallejos/BASiCS/wiki."
  
  msg <- strwrap(msg, exdent=4, indent=4)
  
  packageStartupMessage(paste(msg, collapse="\n"), appendLF=TRUE)
}