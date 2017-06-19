# ==========================================================================
# package initialization
# ==========================================================================
.onAttach = function(libname, pkgname) {
  msg = "BASiCS (v0.5.7 onwards) has been extended to cover differential expression analyses (mean and over-dispersion), as described in http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0930-3. New vignette file is in preparation. In the meantime, please refer to the supplementary files provided in the link for the code required to use the new functionalities. If you have any doubts, please contact cnvallej@uc.cl. "
  msg = strwrap(msg, exdent=4, indent=4)
  packageStartupMessage(paste(msg, collapse="\n"), appendLF=TRUE)
}

# Add message describing which functions changed
# And message describing which functions/classes were removed
