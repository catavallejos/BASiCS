# ==============================================================
# package initialization
# ==============================================================
.onAttach = function(libname, pkgname) {
  msg = "Welcome to 'BASiCS'. If you used 'BASiCS' before its release in
Bioconductor, please visit: https://github.com/catavallejos/BASiCS/wiki."

  msg <- strwrap(msg, exdent=4, indent=4)
  
  packageStartupMessage(paste(msg, collapse="\n"), appendLF=TRUE)
}