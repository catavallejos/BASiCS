#ifndef UTIL_H
#define UTIL_H

#include <RcppArmadillo.h>

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;

/* Auxiliary function that converts Rcpp::NumericVector elements into arma::vec elements */
inline arma::vec as_arma(
  NumericVector& x) /* Vector to be converted into an arma::vec class */
{
    return arma::vec(x.begin(), x.length(), false);
}

/* Auxiliary function that converts Rcpp::NumericMatrix elements into arma::mac elements */
inline arma::mat as_arma(NumericMatrix& x) /* Matrix to be converted into an arma::mat class */
{
    return arma::mat(x.begin(), x.nrow(), x.ncol(), false);
}

#endif
