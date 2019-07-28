#include "MCMCcpp.h"

#ifndef UTILS_H
#define UTILS_H

/* Declarations for general utility functions
 * Stored in: general_utils.cpp
 */
#define ZTOL sqrt(DOUBLE_EPS)

double gig_y_gfn(double y, double m, double beta, double lambda);
double rinvgauss(const double mu, const double lambda);
double zeroin_gig(double ax,double bx,
                  double f(double x, double m, double beta, double lambda),
                  double tol,double m,double beta,double lambda);
NumericVector Rgig(const int n, 
                   const double lambda, 
                   const double chi, 
                   const double psi);
double RgigDouble(const double lambda, const double chi, const double psi);
void StartSampler(int const& N);
void EndBurn();
void CurrentIter(int const& i, int const& N);
void EndSampler(int const& N);
void ReportAR(arma::vec const& AR, 
              std::string const& Param);
arma::vec as_arma(NumericVector& x);
arma::mat as_arma(NumericMatrix& x);
arma::mat lgamma_cpp(arma::mat const& x);
arma::vec lgamma_cpp_vec(arma::vec const& x);
double log_sum_exp_cpp(arma::vec const& x_arma);
arma::vec DegubInd(arma::vec ind,
                   int const q,
                   arma::vec const& u, 
                   arma::vec const& log_aux,
                   arma::vec const& y,
                   double const& threshold,
                   std::string const& param);

// [[Rcpp::export]]
arma::vec Hidden_rDirichlet(arma::vec alpha);
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);

/* Declarations for functions used by the main MCMC sampler
 * Stored in: utils_MCMCcpp.cpp
 */
// [[Rcpp::export]]
arma::mat Hidden_muUpdate(
    arma::vec const& mu0, 
    arma::vec const& prop_var, 
    arma::mat const& Counts, 
    arma::vec const& invdelta, 
    arma::vec const& phinu, 
    arma::vec const& sum_bycell_bio, 
    double const& s2_mu,
    int const& q0,
    int const& n,
    arma::vec & mu1,
    arma::vec & u, 
    arma::vec & ind,
    double const& mintol);

arma::mat deltaUpdate(
    arma::vec const& delta0, 
    arma::vec const& prop_var,  
    arma::mat const& Counts, 
    arma::vec const& mu, 
    arma::vec const& phinu, 
    double const& a_delta, 
    double const& b_delta, 
    double const& s2delta,
    double const& prior_delta,
    int const& q0,
    int const& n,
    arma::vec & delta1,
    arma::vec & u, 
    arma::vec & ind,
    double const& mintol);

Rcpp::List phiUpdate(
    arma::vec const& phi0, 
    double const& prop_var, 
    arma::mat const& Counts, 
    arma::vec const& mu, 
    arma::vec const& invdelta, 
    arma::vec const& nu, 
    arma::vec const& aphi, 
    arma::vec const& sum_bygene_bio, 
    int const& q0,
    int const& n,
    arma::vec & phi1);

arma::vec sUpdateBatch(
    arma::vec const& s0, 
    arma::vec const& nu, 
    arma::vec const& thetaBatch, 
    double const& as, 
    double const& bs, 
    arma::mat const& BatchDesign,
    int const& n,
    arma::vec & s1);

arma::mat nuUpdateBatch(
    arma::vec const& nu0, 
    arma::vec const& prop_var, 
    arma::mat const& Counts,
    double const& SumSpikeInput,
    arma::mat const& BatchDesign, 
    arma::vec const& mu, 
    arma::vec const& invdelta, 
    arma::vec const& phi, 
    arma::vec const& s,
    arma::vec const& thetaBatch, 
    arma::vec const& sum_bygene_all, 
    int const& q0,
    int const& n,
    arma::vec & nu1,
    arma::vec & u,
    arma::vec & ind,
    double const& mintol);

arma::mat thetaUpdateBatch(
    arma::vec const& theta0, 
    arma::vec const& prop_var, 
    arma::mat const& BatchDesign,
    arma::vec const& BatchSizes,
    arma::vec const& s, 
    arma::vec const& nu, 
    double const& a_theta, 
    double const& b_theta, 
    int const& n,
    int const& nBatch,
    double const& mintol);
  
  /* Declarations for functions used by the MCMC sampler for the regression case
   * Stored in: utils_MCMCcppReg.cpp
   */
  
arma::mat designMatrix(
    int const& k, 
    arma::vec const& mu, 
    double const& variance);

arma::mat muUpdateReg(
    arma::vec const& mu0, 
    arma::vec const& prop_var, 
    arma::mat const& Counts, 
    arma::vec const& delta, 
    arma::vec const& phinu, 
    arma::vec const& sum_bycell_bio, 
    double const& s2_mu,
    int const& q0,
    int const& n,
    arma::vec & mu1,
    arma::vec & u, 
    arma::vec & ind,
    int const& k,
    arma::vec const& lambda,
    arma::vec const& beta,
    arma::mat const& X,
    double const& sigma2,
    double variance,
    double const& mintol);

arma::mat deltaUpdateReg(
    arma::vec const& delta0, 
    arma::vec const& prop_var,  
    arma::mat const& Counts, 
    arma::vec const& mu, 
    arma::vec const& phinu, 
    int const& q0,
    int const& n,
    arma::vec & delta1,
    arma::vec & u, 
    arma::vec & ind,
    arma::vec const& lambda,
    arma::mat const& X,
    double const& sigma2,
    arma::vec const& beta,
    double const& mintol);

arma::vec betaUpdateReg(double const& sigma2,
                        arma::mat const& VAux,
                        arma::vec const& mAux);

double sigma2UpdateReg(arma::vec const& delta,
                       arma::vec const& beta,
                       arma::vec const& lambda, 
                       arma::mat const& V1,
                       double const& mInvVm0,
                       arma::vec const& m,
                       double const& sigma2_a0,
                       double const& sigma2_b0,
                       int const& q0);

arma::vec lambdaUpdateReg(arma::vec const& delta,
                          arma::mat const& X,
                          arma::vec const& beta,
                          double const& sigma2, 
                          double const& eta, 
                          int const& q0,
                          arma::vec lambda1);

/* Declarations for functions used by the MCMC sampler for the non-spikes case
 * Stored in: utils_MCMCcppNoSpikes.cpp
 */

arma::mat muUpdateNoSpikes(
    arma::vec const& mu0, 
    arma::vec const& prop_var, 
    arma::mat const& Counts,  
    arma::vec const& invdelta, 
    arma::vec const& nu, 
    arma::vec const& sum_bycell_all, 
    double const& s2_mu,
    int const& q0,
    int const& n,
    arma::vec & mu1,
    arma::vec & u,
    arma::vec & ind,
    double const& Constrain,
    int const& RefGene,
    arma::uvec const& ConstrainGene,
    arma::uvec const& NotConstrainGene,
    int const& ConstrainType,
    double const& mintol);

arma::mat nuUpdateBatchNoSpikes(
    arma::vec const& nu0, 
    arma::vec const& prop_var, 
    arma::mat const& Counts,
    arma::mat const& BatchDesign, 
    arma::vec const& mu, 
    arma::vec const& invdelta, 
    arma::vec const& s, 
    arma::vec const& thetaBatch, 
    arma::vec const& sum_bygene_all, 
    int const& q0,
    int const& n,
    arma::vec & nu1,
    arma::vec & u,
    arma::vec & ind,
    double const& mintol);

/* Declarations for functions used by the MCMC sampler for the regression and 
 * non-spikes case
 * Stored in: utils_MCMCcppRegNoSpikes.cpp
 */

arma::mat muUpdateRegNoSpikes(
    arma::vec const& mu0, 
    arma::vec const& prop_var, 
    arma::mat const& Counts,  
    arma::vec const& delta,
    arma::vec const& invdelta, 
    arma::vec const& nu, 
    arma::vec const& sum_bycell_all, 
    double const& s2_mu,
    int const& q0,
    int const& n,
    arma::vec & mu1,
    arma::vec & u,
    arma::vec & ind,
    double const& Constrain, 
    int const& RefGene,
    arma::uvec const& ConstrainGene,
    arma::uvec const& NotConstrainGene,
    int const& ConstrainType,
    int const& k, 
    arma::vec const& lambda,
    arma::vec const& beta,
    arma::mat const& X,
    double const& sigma2,
    double variance,
    double const& mintol);

arma::mat deltaUpdateRegNoSpikes(
    arma::vec const& delta0,
    arma::vec const& prop_var,  
    arma::mat const& Counts, 
    arma::vec const& mu, 
    arma::vec const& nu, 
    int const& q0,
    int const& n,
    arma::vec & delta1,
    arma::vec & u, 
    arma::vec & ind,
    arma::vec const& lambda,
    arma::mat const& X,
    double const& sigma2,
    arma::vec const& beta,
    double const& mintol);
  
#endif