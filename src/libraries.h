#ifndef MCMCCPP_H
#define MCMCCPP_H

#include <math.h>
#include <assert.h>
#include <R.h>
#include <Rmath.h>
#include <cmath>
//File: matprod_arma.cpp
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]


using namespace Rcpp;

#endif
