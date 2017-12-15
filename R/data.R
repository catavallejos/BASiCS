#' Extract from the chain obtained in Vallejos et al (2016): pool-and-split samples
#'
#' Small extract (100 MCMC iterations, 500 randomly selected genes) from the 
#' chain obtained in Vallejos et al (2016), related to pool-and-split samples 
#' (this corresponds to the RNA 2i samples in Grun et al, 2014).
#' 
#' @format An object of class \code{\linkS4class{BASiCS_Chain}} 
#' containing 100 MCMC iterations.
#'
#' @references 
#' 
#' Vallejos, Richardson and Marioni (2016). Genome Biology.
#' 
#' Grun, Kester and van Oudenaarden (2014). Nature Methods.
#' 
"ChainRNA"

#' Extract from the chain obtained in Vallejos et al (2016): single-cell samples
#'
#' Small extract (100 MCMC iterations, 500 randomly selected genes) from the 
#' chain obtained in Vallejos et al (2016), related to single-cell samples 
#' (this corresponds to the SC 2i samples in Grun et al, 2014). 
#' 
#' @format An object of class \code{\linkS4class{BASiCS_Chain}} 
#' containing 100 MCMC iterations.
#'
#' @references 
#' 
#' Vallejos, Richardson and Marioni (2016). Genome Biology.
#' 
#' Grun, Kester and van Oudenaarden (2014). Nature Methods. 
#' 
"ChainSC"

#' Extract from the chain obtained from the regression model: pool-and-split samples
#'
#' Small extract (10 MCMC iterations, 1000 randomly selected genes) from the 
#' chain obtained by the regression model, related to pool-and-split samples 
#' (this corresponds to the RNA 2i samples in Grun et al, 2014). 
#' 
#' @format An object of class \code{\linkS4class{BASiCS_Chain}} 
#' containing 10 MCMC iterations.
#'
#' @references 
#' 
#' Grun, Kester and van Oudenaarden (2014). Nature Methods. 
#' 
"ChainRNAReg"

#' Extract from the chain obtained from the regression model: single-cell samples
#'
#' Small extract (10 MCMC iterations, 1000 randomly selected genes) from the 
#' chain obtained by the regression model, related to single-cell samples 
#' (this corresponds to the SC 2i samples in Grun et al, 2014). 
#' 
#' @format An object of class \code{\linkS4class{BASiCS_Chain}} 
#' containing 10 MCMC iterations.
#'
#' @references 
#' 
#' Grun, Kester and van Oudenaarden (2014). Nature Methods. 
#' 
"ChainSCReg"
