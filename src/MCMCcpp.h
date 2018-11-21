/* C++ implementation of BASiCS 
 * Author: Catalina A. Vallejos (cnvallej@uc.cl) & Nils Eling
 */

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

using namespace Rcpp;

// Main MCMC sampler
Rcpp::List HiddenBASiCS_MCMCcpp(
    int N, 
    int Thin, 
    int Burn,  
    NumericMatrix Counts,  
    NumericMatrix BatchDesign, 
    NumericVector muSpikes, 
    NumericVector mu0, 
    NumericVector delta0, 
    NumericVector phi0, 
    NumericVector s0,
    NumericVector nu0, 
    NumericVector theta0, 
    double s2mu,
    double adelta, 
    double bdelta, 
    double s2delta,
    double prior_delta,
    NumericVector aphi, 
    double as, 
    double bs,   
    double atheta, 
    double btheta, 
    double ar, 
    NumericVector LSmu0, 
    NumericVector LSdelta0, 
    double LSphi0, 
    NumericVector LSnu0, 
    NumericVector LStheta0, 
    NumericVector sumByCellAll, 
    NumericVector sumByCellBio,
    NumericVector sumByGeneAll, 
    NumericVector sumByGeneBio,
    int StoreAdapt, 
    int EndAdapt,
    int PrintProgress,
    double geneExponent,
    double cellExponent);

// MCMC sampler for regression case
Rcpp::List HiddenBASiCS_MCMCcppReg(
    int N, 
    int Thin, 
    int Burn,  
    NumericMatrix Counts,  
    NumericMatrix BatchDesign, 
    NumericVector muSpikes, 
    NumericVector mu0, 
    NumericVector delta0, 
    NumericVector phi0, 
    NumericVector s0,
    NumericVector nu0, 
    NumericVector theta0, 
    double s2mu,
    NumericVector aphi, 
    double as, 
    double bs,   
    double atheta, 
    double btheta, 
    double ar, 
    NumericVector LSmu0, 
    NumericVector LSdelta0, 
    double LSphi0, 
    NumericVector LSnu0, 
    NumericVector LStheta0, 
    NumericVector sumByCellAll, 
    NumericVector sumByCellBio,
    NumericVector sumByGeneAll, 
    NumericVector sumByGeneBio,
    int StoreAdapt, 
    int EndAdapt,
    int PrintProgress,
    int k,
    NumericVector m0, 
    NumericMatrix V0, 
    double sigma2_a0, 
    double sigma2_b0,
    NumericVector beta0, 
    double sigma20, 
    double eta0, 
    NumericVector lambda0, 
    double const& variance);

// MCMC sampler for the non-spike case
Rcpp::List HiddenBASiCS_MCMCcppNoSpikes(
    int N, 
    int Thin, 
    int Burn,  
    NumericMatrix Counts, 
    NumericMatrix BatchDesign, 
    NumericVector mu0, 
    NumericVector delta0,   
    NumericVector s0, 
    NumericVector nu0,    
    NumericVector theta0,
    double s2mu, 
    double adelta, 
    double bdelta, 
    double s2delta,
    double prior_delta,
    double as,
    double bs,
    double atheta, 
    double btheta, 
    double ar, 
    NumericVector LSmu0, 
    NumericVector LSdelta0, 
    NumericVector LSnu0, 
    NumericVector LStheta0, 
    NumericVector sumByCellAll, 
    NumericVector sumByGeneAll, 
    int StoreAdapt, 
    int EndAdapt,
    int PrintProgress,
    double Constrain,
    NumericVector Index,
    int RefGene,
    NumericVector RefGenes,
    NumericVector ConstrainGene,
    NumericVector NotConstrainGene,
    int ConstrainType,
    int StochasticRef);

// MCMC sampler for regression and non-spikes case
Rcpp::List HiddenBASiCS_MCMCcppRegNoSpikes(
    int N, 
    int Thin, 
    int Burn,  
    NumericMatrix Counts,  
    NumericMatrix BatchDesign, 
    NumericVector mu0, 
    NumericVector delta0, 
    NumericVector s0, 
    NumericVector nu0, 
    NumericVector theta0, 
    double s2mu,
    double as,
    double bs,
    double atheta, 
    double btheta, 
    double ar, 
    NumericVector LSmu0, 
    NumericVector LSdelta0, 
    NumericVector LSnu0, 
    NumericVector LStheta0, 
    NumericVector sumByCellAll, 
    NumericVector sumByGeneAll, 
    int StoreAdapt, 
    int EndAdapt,
    int PrintProgress,
    int k,
    NumericVector m0, 
    NumericMatrix V0, 
    double sigma2_a0, 
    double sigma2_b0,
    NumericVector beta0, 
    double sigma20, 
    double eta0, 
    NumericVector lambda0, 
    double const& variance,
    double Constrain,
    NumericVector Index,
    int RefGene,
    NumericVector RefGenes,
    NumericVector ConstrainGene,
    NumericVector NotConstrainGene,
    int ConstrainType,
    int StochasticRef);

// Function to cumpute denoised rates
arma::mat HiddenBASiCS_DenoisedRates(
    NumericMatrix CountsBio, 
    NumericMatrix Mu,
    NumericMatrix TransInvDelta,
    NumericMatrix PhiNu, 
    int N,
    int q0,
    int n);

#endif
